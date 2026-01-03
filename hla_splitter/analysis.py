import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.manifold import TSNE
import scipy.spatial as spt
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import seaborn as sns
import warnings
import logging
import os

log = logging.getLogger(__name__)

def prepare_hla_inputs(hla_list_path, out_dir):
    """
    Processes the input HLA list CSV to extract unique genotypes.
    Equivalent to HLA_input.py.
    """
    log.info(f"Processing HLA list: {hla_list_path}")
    try:
        sampleInfo = pd.read_csv(hla_list_path, header=0, index_col=0)
        HLA_genotype = np.unique(np.array(sampleInfo))
        out_file = os.path.join(out_dir, "sorted_genotypes.txt")
        pd.DataFrame(HLA_genotype).to_csv(out_file, header=False, index=False)
        log.info(f"Unique HLA genotypes saved to {out_file}")
        return out_file
    except FileNotFoundError:
        log.error(f"HLA list file not found: {hla_list_path}")
        raise
    except Exception as e:
        log.error(f"Error processing HLA list: {e}")
        raise

def summarize_counts(out_dir):
    """
    Summarizes bustools counts per HLA gene based on sorted genotypes.
    Equivalent to Summarize.py.
    """
    log.info("Summarizing bustools counts per HLA type...")
    try:
        filter_counts_dir = os.path.join(out_dir, "sc_output", "filter_counts")
        barcodes_file = os.path.join(filter_counts_dir, "output.barcodes.txt")
        genes_file = os.path.join(filter_counts_dir, "output.genes.txt")
        mtx_file = os.path.join(filter_counts_dir, "output.mtx")
        hla_genotypes_file = os.path.join(out_dir, "sorted_genotypes.txt")
        hla_type_mtx_file = os.path.join(filter_counts_dir, "HLA_type.mtx")

        barcodes = pd.read_csv(barcodes_file, header=None).iloc[:, 0]
        genes = pd.read_csv(genes_file, header=None).iloc[:, 0]
        # Prefix genes with 'HLA-' to match expected format if needed
        # This depends on how t2g was created - adjust if necessary
        # genes = "HLA-" + genes

        HLA_labels = pd.read_csv(hla_genotypes_file, header=None).iloc[:, 0]
        # Prefix HLA labels with 'HLA-' if genes have it
        # HLA_labels = "HLA-" + HLA_labels

        log.info("Reading MTX file...")
        with open(mtx_file, "r") as fileHandler:
            header = [next(fileHandler) for _ in range(3)] # Read header lines
            # Efficiently read data using pandas or iterate line by line for large files
            data = pd.read_csv(fileHandler, sep=' ', skiprows=0, header=None, comment='%',
                               names=['cell_num', 'gene_num', 'count']) # Note: MTX is cell, gene,count
        log.info("Finished reading MTX file.")

        data['cell_num'] = pd.to_numeric(data['cell_num'])
        data['gene_num'] = pd.to_numeric(data['gene_num'])
        data['count'] = pd.to_numeric(data['count'])

        # Map gene numbers (1-based in MTX) to gene names (0-based pandas index)
        gene_map = {i + 1: name for i, name in enumerate(genes)}
        data['gene_name'] = data['gene_num'].map(gene_map)

        data_sum = []
        log.info("Aggregating counts per HLA type...")
        for i, hla_label in enumerate(HLA_labels):
            # Find genes matching the current HLA label (might need flexible matching)
            # Assuming exact match for now, adjust if partial matching needed
            matching_genes = genes[genes.astype(str).str.startswith(hla_label)]
            if matching_genes.empty:
                log.warning(f"No matching gene found in counts for HLA type: {hla_label}")
                continue

            # Get the original 1-based indices for matching genes
            gene_indices = matching_genes.index + 1
            data_tmp = data[data['gene_num'].isin(gene_indices)]

            # Group by cell and sum counts for the current HLA type
            hla_sum_per_cell = data_tmp.groupby('cell_num')['count'].sum().round(2)

            for cell_idx, total_count in hla_sum_per_cell.items():
                # Store cell index (1-based), HLA type index (1-based), and summed count
                data_sum.append([cell_idx, i + 1, total_count])

            log.info(f"Processed HLA type: {hla_label} ({len(hla_sum_per_cell)} cells)")

        data_sum.sort() # Sort primarily by cell index

        log.info(f"Writing summarized HLA matrix to {hla_type_mtx_file}")
        with open(hla_type_mtx_file, "w") as f:
            f.write(header[0]) # Banner
            f.write(header[1]) # Comment
            # Write updated header: features, cells, non-zero entries
            f.write(f"{len(barcodes)} {len(HLA_labels)} {len(data_sum)}\n") # Correct order: features cells entries
            for entry in data_sum:
                 # Output: cell_idx count hla_idx(feature) 
                 f.write(f"{entry[0]} {entry[1]} {entry[2]}\n")

        log.info("Summarizing counts finished.")
        return hla_type_mtx_file
    except FileNotFoundError as e:
        log.error(f"Required file not found during summarization: {e}")
        raise
    except Exception as e:
        log.error(f"Error during count summarization: {e}")
        raise


def run_demultiplex(hla_list_path, out_dir):
    """
    Performs cell demultiplexing based on summarized HLA counts and HLA list.
    Equivalent to Demltiplex.py.
    """
    log.info("Starting HLA-Splitter...")
    try:
        filter_counts_dir = os.path.join(out_dir, "sc_output", "filter_counts")
        hla_genotypes_file = os.path.join(out_dir, "sorted_genotypes.txt")
        barcodes_file = os.path.join(filter_counts_dir, "output.barcodes.txt")
        hla_type_mtx_file = os.path.join(filter_counts_dir, "HLA_type.mtx")
        output_csv = os.path.join(out_dir, "HLA_demultiplex.csv")

        log.info("Loading input data...")
        sampleInfo = pd.read_csv(hla_list_path, header=0, index_col=0)
        namelist = sampleInfo.columns # 获取样本名称列表
        HLA_list_per_sample = [pd.unique(sampleInfo[col]) for col in namelist]
        HLA_labels = pd.read_csv(hla_genotypes_file, header=None).iloc[:, 0]
        barcodes = pd.read_csv(barcodes_file, header=None).iloc[:, 0]

        # --- 1. 数据读取与准备 ---
        log.info(f"Reading summarized HLA matrix: {hla_type_mtx_file}")
        try:
            adata = sc.read_mtx(hla_type_mtx_file)
            adata.var_names = HLA_labels  # Features are HLA types
            adata.obs_names = barcodes    # Observations are cells
            HLA_matrix = adata.to_df()    # Convert to dense DataFrame (obs x vars)
            log.info(f"Matrix shape: {HLA_matrix.shape}")
        except Exception as e:
            log.error(f"Error reading HLA matrix or assigning labels: {e}")

        # --- 2. 归一化与缩放 ---
        log.info("Performing total count normalization and log1p transformation...")
    
        # 转换为 AnnData 对象进行归一化，因为 scanpy 的功能很方便
        # 确保 AnnData 的 X 是 float 类型，以避免归一化问题
        adata_hla = sc.AnnData(HLA_matrix.astype(float)) 
        adata_hla.obs_names = HLA_matrix.index
        adata_hla.var_names = HLA_matrix.columns

        # 步骤 2.1: 总计数归一化 (Normalize total counts to 100 per cell)
        sc.pp.normalize_total(adata_hla, target_sum=1e2) 

        # 步骤 2.2: 对数转换 (log1p transformation: log(x+1))
        sc.pp.log1p(adata_hla)

        # 步骤 2.3: 特征缩放 (Scale each feature/HLA type to unit variance, usually after log transform)
        log.info("Scaling features (HLA types) to unit variance...")
    
        # 手动实现类似原始的对 positive 值缩放，但更规范
        HLA_matrix_norm_scaled = adata_hla.to_df().copy()
        scaler = StandardScaler(with_mean=False, with_std=True)
        warnings.filterwarnings("ignore", category=UserWarning, module='sklearn.preprocessing._data')
        for hla_type in HLA_matrix_norm_scaled.columns:
            exp = HLA_matrix_norm_scaled[hla_type]
            exp_positive_indices = exp[exp > 0].index
            if not exp_positive_indices.empty:
                scaled_values = scaler.fit_transform(exp[exp_positive_indices].values.reshape(-1, 1)).flatten()
                HLA_matrix_norm_scaled.loc[exp_positive_indices,hla_type] = scaled_values
            else:
                log.debug(f"HLA type {hla_type} has no positive expression after log1p, skipping scaling.")
    
        HLA_matrix_final = HLA_matrix_norm_scaled
    
    
        # --- 3. 计算样本得分 ---
        log.info("Calculating demultiplexing scores...")
        score = pd.DataFrame(index=HLA_matrix_final.index)
        HLA_types_major = ["A", "B", "C", "DRB1", "DQB1", "DPB1"] # Major HLA loci

        for i, sample_name in enumerate(namelist):
            score_tmp = pd.Series(0.0, index=HLA_matrix_final.index)
            sample_alleles = HLA_list_per_sample[i]
            log.debug(f"Processing sample: {sample_name} with alleles: {sample_alleles}")
        
            for hla_locus in HLA_types_major:
                locus_alleles = [allele for allele in sample_alleles if allele.startswith(hla_locus)]
                # 纯合子加权因子，如果只有一个等位基因，则视为纯合
                homo_factor = 1.8 if len(locus_alleles) == 1 else 1
            
                for allele in locus_alleles:
                    if allele in HLA_matrix_final.columns:
                        tmp = HLA_matrix_final[allele]
                        # 只加有表达的等位基因，mean() != 0 可以简化为直接判断是否有值，因为tmp可能是稀疏的
                        # 或者如果tmp的sum() > 0
                        if tmp.sum() > 0: # 检查是否有非零表达
                            score_tmp += tmp * homo_factor
                        else:
                            log.debug(f"Allele {allele} has zero sum expression across cells for scoring.")
                    else:
                        log.debug(f"Allele {allele} not found in HLA matrix columns.")
            score[f"{sample_name}"] = score_tmp
        
        # 初始预测基于最大得分
        log.info("Assigning initial predictions...")
        sample_predict_series = score.idxmax(axis=1)
        # 处理所有得分都为零的细胞
        sample_predict_series[score.sum(axis=1) == 0] = "Undefined"
 
        # --- 4. t-SNE 嵌入 ---
        log.info("Calculating t-SNE embedding...")
        if HLA_matrix_final.shape[0] < 2 or HLA_matrix_final.shape[1] < 2:
            log.warning("Matrix too small for t-SNE, skipping embedding and adjustment. Returning initial predictions.")
            tsne_df = pd.DataFrame(index=HLA_matrix_final.index)
            tsne_df["sample_predict"] = sample_predict_series
            tsne_df["adjust_predict"] = sample_predict_series
            certainty = pd.Series(1.0, index=HLA_matrix_final.index) # 默认高置信度
        else:
            # 调整 perplexity，确保其小于细胞总数-1
            # perplexity 至少为 1，通常在 5-50 之间
            perplexity_val = min(max(5, HLA_matrix_final.shape[0] - 1), 50) # 确保 perplexity 在合理范围内
            log.info(f"Using t-SNE perplexity: {perplexity_val}")
        
            # 检查 NaN/Inf 值，t-SNE不能处理它们
            if HLA_matrix_final.isnull().any().any() or np.isinf(HLA_matrix_final).any().any():
                log.warning("HLA matrix contains NaN or Inf values. Attempting to fill with 0 for t-SNE.")
                # 填充 NaN/Inf，或者考虑更复杂的数据清洗
                HLA_matrix_for_tsne = HLA_matrix_final.fillna(0).replace([np.inf, -np.inf], 0)
            else:
                HLA_matrix_for_tsne = HLA_matrix_final

            tsne = TSNE(n_components=2, init="pca", random_state=42, perplexity=perplexity_val, n_jobs=-1) # n_jobs=-1 使用所有可用CPU
            tsne_embedding = tsne.fit_transform(HLA_matrix_for_tsne.values)
            tsne_df = pd.DataFrame(tsne_embedding, index=HLA_matrix_for_tsne.index, columns=['TSNE1', 'TSNE2'])
            tsne_df["sample_predict"] = sample_predict_series

        # --- 5. 基于邻域的预测调整 ---
        log.info("Adjusting predictions based on neighbors...")
        n_neighbor = min(100, HLA_matrix_final.shape[0] - 1) # 确保 n_neighbor 不超过实际细胞数-1
        if n_neighbor > 1: # 至少需要2个点才能计算邻居
            point_set = tsne_df[['TSNE1', 'TSNE2']].values
            log.info("Building KDTree...")
            kt = spt.KDTree(point_set)
            log.info(f"Querying {n_neighbor} neighbors...")
            _, neighbour_indices = kt.query(point_set, k=n_neighbor + 1) # 查询 k+1 个，第一个是自身

            log.info("Adjusting labels based on neighbors...")
            adjusted_predict = tsne_df["sample_predict"].copy()
            certainty_list = []

            for i in range(len(tsne_df)):
                valid_neighbor_indices = [idx for idx in neighbour_indices[i] if idx != i]
                if not valid_neighbor_indices: # 如果没有有效邻居（例如，细胞总数太少）
                    certainty_list.append(0.0)
                    continue
    
                neighbor_labels = tsne_df["sample_predict"].iloc[valid_neighbor_indices]
                neighbor_counts = neighbor_labels.value_counts()
    
                if neighbor_counts.empty:
                    certainty_list.append(0.0)
                    continue
    
                dominate_type = neighbor_counts.index[0]
                dominate_count = neighbor_counts.iloc[0]
                cert = dominate_count / (len(valid_neighbor_indices)) # 计算置信度，基于实际邻居数
                certainty_list.append(cert)
    
                current_label = tsne_df["sample_predict"].iloc[i]
    
                if dominate_type != current_label:
                    # 如果邻居主导类型的置信度足够高，且当前细胞的原始分数不高，则调整
                    if cert > 0.5: # 可以根据需要调整这个阈值
                        adjusted_predict.iloc[i] = dominate_type
                    else:
                        adjusted_predict.iloc[i] = "Undefined"
                
            tsne_df["adjust_predict"] = adjusted_predict
            certainty = pd.Series(certainty_list, index=tsne_df.index)
        else:
            log.warning(f"Not enough cells ({HLA_matrix_final.shape[0]}) for neighbor adjustment with {n_neighbor} neighbors. Skipping neighbor adjustment.")
            tsne_df["adjust_predict"] = tsne_df["sample_predict"]
            certainty = pd.Series(1.0, index=tsne_df.index)

        # --- 6. 绘图：t-SNE 散点图 ---
        # 确保 category_list 包含所有可能的类别，包括 "Undefined"
        # namelist 必须是 list 或可迭代对象，不是 pd.Index
        if isinstance(namelist, pd.Index):
            category_list_plot = namelist.tolist()
        else:
            category_list_plot = list(namelist)
        if "Undefined" not in category_list_plot:
            category_list_plot.append("Undefined")

        tsne_df['adjust_predict'] = pd.Categorical(
            tsne_df['adjust_predict'],
            categories=category_list_plot,
            ordered=True
        )

        actual_colored_samples = tsne_df[tsne_df['adjust_predict'] != "Undefined"]['adjust_predict'].unique()
        n_colors_for_legend = len(actual_colored_samples)

        if n_colors_for_legend > 0 or "Undefined" in tsne_df['adjust_predict'].unique(): 
            # 增加图形宽度，为图例留出更多空间
            plt.figure(figsize=(10, 8)) # 10x8 inches 适合大多数情况
        
            # --- Create a custom color palette ---
            base_palette = sns.color_palette("deep", n_colors=max(1, len(namelist)))
            # Create a dictionary to map sample names to colors
            custom_palette = {}
            for i, sample_name in enumerate(namelist): # Sort to ensure consistent color assignment
                if i < len(base_palette): # Assign colors from the base palette
                    custom_palette[sample_name] = base_palette[i]
                else: # Fallback if more samples than default palette colors (should ideally not happen with "deep")
                    custom_palette[sample_name] = np.random.rand(3,) # Assign a random color if palette runs out
    
            # Assign gray to "Undefined"
            custom_palette["Undefined"] = "gray"
    
            # 绘制散点图
            sns.scatterplot(x='TSNE1', y='TSNE2', hue='adjust_predict', data=tsne_df[tsne_df['adjust_predict'] != "Undefined"], # 使用所有数据，包括 Undefined
                            palette=custom_palette, s=3, linewidth=0)
            plt.title('HLA t-SNE with Sample labels')

            # 绘制质心标签
            # 仅对非 "Undefined" 的样本计算和绘制质心
            centroids_adjusted = tsne_df[tsne_df['adjust_predict']!="Undefined"].groupby('adjust_predict')[['TSNE1', 'TSNE2']].mean().reset_index()
            for index, row in centroids_adjusted.iterrows():
                plt.text(row['TSNE1'], row['TSNE2'], row['adjust_predict'], fontsize=10, color='black',
                         ha='center', va='center') # 居中对齐标签
            
            plt.xlabel('tSNE1')
            plt.ylabel('tSNE2')
        
            # 动态计算图例列数，确保至少为1列，并限制最大列数（例如，不超过5列）
            # 调整分母，使其更合理。如果图例项很多，可以考虑更大的分母，让每列更短。
            if n_colors_for_legend > 0:
                legend_ncol = max(1, round(n_colors_for_legend / 20)) # 每列大约20个项
                legend_ncol = min(legend_ncol, 3) # 限制最大列数，防止图例过宽
            else:
                legend_ncol = 1 # 如果没有可着色样本，只有 Undefined 或少数几个，保持1列

            plt.legend(title='Adjusted_Prediction', bbox_to_anchor=(1.02, 1), loc='upper left', ncol=legend_ncol, borderaxespad=0.)
            plt.tight_layout()
            plt.savefig(os.path.join(out_dir, "HLA_TSNE_Adjusted_Prediction.pdf"))
            plt.close()
        else:
            log.warning("No valid samples for t-SNE plot after 'Undefined' filtering. Skipping plot.")


        # --- 7. 绘图：Countplot ---
        log.info("Generating countplot of adjusted predictions...")
        plt.figure(figsize=(10, 6))
        # 确保排序，并可以根据需要调整 x 轴标签的顺序
        sns.countplot(data=tsne_df, x='adjust_predict', palette=custom_palette, order=category_list_plot)
        plt.xlabel("Adjusted Prediction")
        plt.ylabel("Cell numbers")
        plt.title("Distribution of Sample Label")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, "Countplot_Adjusted_Prediction.pdf"))
        plt.close()

        # --- 8. 最终输出 ---
        log.info(f"Saving demultiplexing results to {output_csv}")
        predict_result = pd.DataFrame({
            "Cell barcode": tsne_df.index,
            "Prediction": tsne_df["sample_predict"],
            "Adjusted_prediction": tsne_df["adjust_predict"],
            "Certainty": certainty
        })
        predict_result.to_csv(output_csv, index=False)
        log.info("Demultiplexing finished.")
        return output_csv

    except FileNotFoundError as e:
        log.error(f"Required file not found during demultiplexing: {e}")
        raise
    except Exception as e:
        log.error(f"Error during demultiplexing: {e}")
        raise