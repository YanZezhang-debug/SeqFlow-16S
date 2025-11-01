#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
基于图论方法的de novo组装算法
使用重叠图（overlap graph）处理Sanger测序数据
"""

import os
import re
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner
import numpy as np
from collections import defaultdict, deque
import logging
import networkx as nx
from typing import Dict, List, Tuple, Set, Optional

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class GraphBasedAssembler:
    def __init__(self, min_overlap=50, min_quality=20, min_length=100):
        """
        初始化基于图论的组装器
        
        Args:
            min_overlap: 最小重叠长度
            min_quality: 最小质量分数
            min_length: 最小序列长度
        """
        self.min_overlap = min_overlap
        self.min_quality = min_quality
        self.min_length = min_length
        self.sequences = {}  # 存储所有序列
        self.qualities = {}  # 存储质量分数
        self.overlap_graph = nx.DiGraph()  # 重叠图
        self.contigs = []  # 组装的contig
        
    def parse_filename(self, filename: str) -> Optional[Dict]:
        """解析文件名，提取样本信息和测序方向"""
        pattern = r'(\d+)_(\d+)_\(([^)]+)\)_\[([^\]]+)\]\.ab1'
        match = re.match(pattern, filename)
        if match:
            file_num, sample_id, sample_name, direction = match.groups()
            return {
                'file_num': file_num,
                'sample_id': sample_id,
                'sample_name': sample_name,
                'direction': direction,
                'filename': filename
            }
        return None
    
    def read_ab1_files(self, data_dir: str) -> Dict:
        """读取所有ab1文件"""
        logger.info("正在读取ab1文件...")
        
        samples = defaultdict(dict)
        ab1_files = [f for f in os.listdir(data_dir) if f.endswith('.ab1')]
        logger.info(f"找到 {len(ab1_files)} 个ab1文件")
        
        for filename in ab1_files:
            file_info = self.parse_filename(filename)
            if file_info:
                sample_name = file_info['sample_name']
                direction = file_info['direction']
                samples[sample_name][direction] = file_info
                logger.info(f"  已解析: {filename} -> 样本={sample_name}, 方向={direction}")
            else:
                logger.warning(f"  无法解析文件名: {filename}")
        
        logger.info(f"成功解析 {len(samples)} 个样本")
        
        # 显示每个样本的测序方向
        for sample_name, directions in samples.items():
            direction_list = list(directions.keys())
            logger.info(f"  样本 {sample_name}: {direction_list}")
        
        return samples
    
    def read_ab1_file(self, filepath: str) -> Optional[Dict]:
        """读取单个ab1文件，如果不存在则尝试读取对应的fasta文件"""
        try:
            # 首先尝试读取ab1文件
            if os.path.exists(filepath):
                record = SeqIO.read(filepath, "abi")
                sequence = str(record.seq)
                qualities = record.letter_annotations.get('phred_quality', [])
                
                return {
                    'sequence': sequence,
                    'qualities': qualities,
                    'record': record
                }
            
            # 如果ab1文件不存在，尝试读取对应的fasta文件
            fasta_filepath = filepath.replace('.ab1', '.fasta')
            if os.path.exists(fasta_filepath):
                logger.info(f"ab1文件不存在，使用fasta文件: {fasta_filepath}")
                record = SeqIO.read(fasta_filepath, "fasta")
                sequence = str(record.seq)
                # fasta文件没有质量分数，使用默认值
                qualities = [30] * len(sequence)  # 假设质量分数为30
                
                return {
                    'sequence': sequence,
                    'qualities': qualities,
                    'record': record
                }
            
            logger.error(f"文件不存在: {filepath} 或 {fasta_filepath}")
            return None
            
        except Exception as e:
            logger.error(f"读取文件 {filepath} 时出错: {e}")
            return None
    
    def quality_trim(self, sequence: str, qualities: List[int]) -> Tuple[str, List[int]]:
        """基于质量值修剪序列"""
        if not qualities or len(qualities) != len(sequence):
            return sequence, qualities
        
        # 从两端修剪低质量区域
        start = 0
        end = len(sequence)
        
        # 从5'端修剪
        for i, qual in enumerate(qualities):
            if qual >= self.min_quality:
                start = i
                break
        
        # 从3'端修剪
        for i in range(len(qualities) - 1, -1, -1):
            if qualities[i] >= self.min_quality:
                end = i + 1
                break
        
        trimmed_seq = sequence[start:end]
        trimmed_qual = qualities[start:end]
        
        if len(trimmed_seq) < self.min_length:
            logger.warning(f"修剪后序列过短: {len(trimmed_seq)} < {self.min_length}")
            return sequence, qualities
        
        return trimmed_seq, trimmed_qual
    
    def reverse_complement(self, sequence: str) -> str:
        """获取反向互补序列"""
        return str(Seq(sequence).reverse_complement())
    
    def trim_low_quality_ends(self, sequence: str, qualities: List[int], min_quality: int = 20) -> Tuple[str, List[int]]:
        """修剪序列两端的低质量碱基"""
        if not qualities or len(qualities) != len(sequence):
            return sequence, qualities
        
        # 从5'端修剪低质量碱基
        start = 0
        for i, qual in enumerate(qualities):
            if qual >= min_quality:
                start = i
                break
        
        # 从3'端修剪低质量碱基
        end = len(sequence)
        for i in range(len(qualities) - 1, -1, -1):
            if qualities[i] >= min_quality:
                end = i + 1
                break
        
        trimmed_seq = sequence[start:end]
        trimmed_qual = qualities[start:end]
        
        return trimmed_seq, trimmed_qual
    
    def assemble_forward_reverse(self, forward_seq: str, reverse_seq: str, 
                                forward_qual: List[int], reverse_qual: List[int]) -> Optional[Dict]:
        """组装单个样本的正向和反向序列"""
        # 计算重叠
        overlap_len, start1, start2 = self.calculate_overlap(forward_seq, reverse_seq)
        
        if overlap_len < self.min_overlap:
            logger.warning(f"重叠长度不足: {overlap_len} < {self.min_overlap}")
            return None
        
        # 组装序列
        assembled_sequence = ""
        assembled_qualities = []
        
        # 添加重叠前的正向序列（修剪低质量碱基）
        forward_before = forward_seq[:start1]
        forward_qual_before = forward_qual[:start1]
        trimmed_forward_before, trimmed_forward_qual_before = self.trim_low_quality_ends(
            forward_before, forward_qual_before
        )
        assembled_sequence += trimmed_forward_before
        assembled_qualities.extend(trimmed_forward_qual_before)
        
        # 处理重叠区域 - 选择质量更高的碱基
        for i in range(overlap_len):
            pos1 = start1 + i
            pos2 = start2 + i
            
            if pos1 < len(forward_qual) and pos2 < len(reverse_qual):
                f_qual = forward_qual[pos1]
                r_qual = reverse_qual[pos2]
                
                if f_qual >= r_qual:
                    assembled_sequence += forward_seq[pos1]
                    assembled_qualities.append(f_qual)
                else:
                    assembled_sequence += reverse_seq[pos2]
                    assembled_qualities.append(r_qual)
            else:
                # 如果超出质量数组范围，使用正向序列
                assembled_sequence += forward_seq[pos1]
                assembled_qualities.append(forward_qual[pos1] if pos1 < len(forward_qual) else 0)
        
        # 添加重叠后的反向序列（修剪低质量碱基）
        remaining_reverse = reverse_seq[start2 + overlap_len:]
        remaining_reverse_qual = reverse_qual[start2 + overlap_len:] if start2 + overlap_len < len(reverse_qual) else []
        trimmed_remaining_reverse, trimmed_remaining_reverse_qual = self.trim_low_quality_ends(
            remaining_reverse, remaining_reverse_qual
        )
        assembled_sequence += trimmed_remaining_reverse
        assembled_qualities.extend(trimmed_remaining_reverse_qual)
        
        return {
            'sequence': assembled_sequence,
            'qualities': assembled_qualities,
            'overlap_length': overlap_len
        }
    
    def calculate_overlap(self, seq1: str, seq2: str) -> Tuple[int, int, int]:
        """
        计算两个序列之间的重叠
        
        Returns:
            (overlap_length, start_pos_seq1, start_pos_seq2)
        """
        # 使用更简单的方法计算重叠
        # 寻找最长公共子串
        max_overlap = 0
        best_start1 = 0
        best_start2 = 0
        
        # 检查seq1的每个位置作为起始点
        for i in range(len(seq1)):
            for j in range(len(seq2)):
                # 计算从位置i和j开始的最大匹配长度
                overlap_len = 0
                pos1, pos2 = i, j
                
                while (pos1 < len(seq1) and pos2 < len(seq2) and 
                       seq1[pos1] == seq2[pos2]):
                    overlap_len += 1
                    pos1 += 1
                    pos2 += 1
                
                if overlap_len > max_overlap:
                    max_overlap = overlap_len
                    best_start1 = i
                    best_start2 = j
        
        if max_overlap < self.min_overlap:
            return 0, 0, 0
        
        return max_overlap, best_start1, best_start2
    
    def build_overlap_graph(self, sequences: Dict[str, str]) -> nx.DiGraph:
        """构建重叠图"""
        logger.info("正在构建重叠图...")
        
        # 添加节点
        for seq_id, sequence in sequences.items():
            self.overlap_graph.add_node(seq_id, sequence=sequence)
        
        # 计算所有序列对之间的重叠
        seq_ids = list(sequences.keys())
        for i in range(len(seq_ids)):
            for j in range(i + 1, len(seq_ids)):
                seq1_id = seq_ids[i]
                seq2_id = seq_ids[j]
                seq1 = sequences[seq1_id]
                seq2 = sequences[seq2_id]
                
                # 计算正向重叠
                overlap_len, start1, start2 = self.calculate_overlap(seq1, seq2)
                if overlap_len >= self.min_overlap:
                    self.overlap_graph.add_edge(seq1_id, seq2_id, 
                                              weight=overlap_len, 
                                              start1=start1, start2=start2)
                
                # 计算反向重叠（seq2的反向互补与seq1）
                seq2_rc = self.reverse_complement(seq2)
                overlap_len_rc, start1_rc, start2_rc = self.calculate_overlap(seq1, seq2_rc)
                if overlap_len_rc >= self.min_overlap:
                    self.overlap_graph.add_edge(seq1_id, seq2_id, 
                                              weight=overlap_len_rc, 
                                              start1=start1_rc, start2=start2_rc,
                                              reverse=True)
        
        logger.info(f"重叠图构建完成，包含 {self.overlap_graph.number_of_nodes()} 个节点和 {self.overlap_graph.number_of_edges()} 条边")
        return self.overlap_graph
    
    def find_contigs(self) -> List[List[str]]:
        """在重叠图中找到contig"""
        logger.info("正在查找contig...")
        
        # 使用强连通分量算法找到contig
        contigs = []
        visited = set()
        
        for node in self.overlap_graph.nodes():
            if node in visited:
                continue
            
            # 使用DFS找到连通分量
            contig = self._dfs_contig(node, visited)
            if len(contig) > 1:  # 至少包含2个序列
                contigs.append(contig)
        
        logger.info(f"找到 {len(contigs)} 个contig")
        return contigs
    
    def _dfs_contig(self, start_node: str, visited: Set[str]) -> List[str]:
        """使用DFS查找contig"""
        contig = []
        stack = [start_node]
        
        while stack:
            node = stack.pop()
            if node in visited:
                continue
            
            visited.add(node)
            contig.append(node)
            
            # 添加所有邻居
            for neighbor in self.overlap_graph.neighbors(node):
                if neighbor not in visited:
                    stack.append(neighbor)
        
        return contig
    
    def assemble_contig(self, contig_nodes: List[str]) -> str:
        """组装单个contig"""
        if len(contig_nodes) == 1:
            return self.overlap_graph.nodes[contig_nodes[0]]['sequence']
        
        # 使用贪心算法组装contig
        assembled_seq = ""
        current_node = contig_nodes[0]
        assembled_seq = self.overlap_graph.nodes[current_node]['sequence']
        
        remaining_nodes = set(contig_nodes[1:])
        
        while remaining_nodes:
            best_next = None
            best_overlap = 0
            best_start = 0
            
            # 找到最佳的下一个序列
            for neighbor in self.overlap_graph.neighbors(current_node):
                if neighbor in remaining_nodes:
                    edge_data = self.overlap_graph.edges[current_node, neighbor]
                    if edge_data['weight'] > best_overlap:
                        best_overlap = edge_data['weight']
                        best_next = neighbor
                        best_start = edge_data['start2']
            
            if best_next is None:
                break
            
            # 添加下一个序列
            next_seq = self.overlap_graph.nodes[best_next]['sequence']
            assembled_seq += next_seq[best_start:]
            remaining_nodes.remove(best_next)
            current_node = best_next
        
        return assembled_seq
    
    def generate_consensus(self, sequences: List[str], qualities: List[List[int]]) -> str:
        """基于质量分数生成共识序列"""
        if not sequences:
            return ""
        
        if len(sequences) == 1:
            return sequences[0]
        
        # 找到最短序列的长度作为参考
        min_length = min(len(seq) for seq in sequences)
        consensus = ""
        
        for pos in range(min_length):
            # 收集该位置的所有碱基和质量分数
            bases = []
            base_qualities = []
            
            for i, seq in enumerate(sequences):
                if pos < len(seq) and pos < len(qualities[i]):
                    bases.append(seq[pos])
                    base_qualities.append(qualities[i][pos])
            
            if not bases:
                continue
            
            # 选择质量最高的碱基
            best_quality = max(base_qualities)
            best_base = bases[base_qualities.index(best_quality)]
            consensus += best_base
        
        return consensus
    
    def process_samples(self, samples: Dict, data_dir: str = ".") -> List[Dict]:
        """处理所有样本，每个样本的正向和反向序列组装成一个序列"""
        logger.info("开始处理样本...")
        
        assembled_samples = []
        
        for sample_name, sample_data in samples.items():
            if '16SF' not in sample_data or '16SR' not in sample_data:
                logger.warning(f"样本 {sample_name} 缺少正向或反向序列")
                continue
            
            # 读取正向和反向序列
            forward_file = os.path.join(data_dir, sample_data['16SF']['filename'])
            reverse_file = os.path.join(data_dir, sample_data['16SR']['filename'])
            
            forward_data = self.read_ab1_file(forward_file)
            reverse_data = self.read_ab1_file(reverse_file)
            
            if not forward_data or not reverse_data:
                continue
            
            # 质量修剪
            forward_seq, forward_qual = self.quality_trim(
                forward_data['sequence'], forward_data['qualities']
            )
            reverse_seq, reverse_qual = self.quality_trim(
                reverse_data['sequence'], reverse_data['qualities']
            )
            
            # 获取反向互补序列，同时反转质量分数
            reverse_rc = self.reverse_complement(reverse_seq)
            reverse_qual_rc = list(reversed(reverse_qual))
            
            logger.info(f"样本 {sample_name}: 正向长度 {len(forward_seq)}, 反向长度 {len(reverse_rc)}")
            
            # 组装这个样本的正向和反向序列
            assembled_seq = self.assemble_forward_reverse(forward_seq, reverse_rc, forward_qual, reverse_qual_rc)
            
            if assembled_seq:
                # 评估组装质量
                quality_metrics = self.evaluate_assembly_quality(
                    assembled_seq['sequence'],
                    assembled_seq['qualities'],
                    assembled_seq['overlap_length'],
                    len(forward_seq),
                    len(reverse_rc)
                )
                
                assembled_samples.append({
                    'sample_name': sample_name,
                    'assembled_sequence': assembled_seq['sequence'],
                    'qualities': assembled_seq['qualities'],
                    'length': len(assembled_seq['sequence']),
                    'overlap_length': assembled_seq['overlap_length'],
                    'forward_length': len(forward_seq),
                    'reverse_length': len(reverse_rc),
                    'quality_score': quality_metrics['overall_quality_score'],
                    'quality_grade': quality_metrics['quality_grade'],
                    'gc_content': quality_metrics['gc_content'],
                    'low_quality_ratio': quality_metrics['low_quality_ratio'],
                    'overlap_ratio': quality_metrics['overlap_ratio'],
                    'average_quality': quality_metrics['average_quality']
                })
                logger.info(f"样本 {sample_name} 组装完成，最终长度: {len(assembled_seq['sequence'])}, 质量评分: {quality_metrics['overall_quality_score']}/100 ({quality_metrics['quality_grade']}级)")
            else:
                logger.warning(f"样本 {sample_name} 组装失败")
        
        self.assembled_samples = assembled_samples
        return assembled_samples
    
    def assemble_all_sequences(self) -> List[Dict]:
        """组装所有序列"""
        # 构建重叠图
        self.build_overlap_graph(self.sequences)
        
        # 找到contig
        contigs = self.find_contigs()
        
        # 组装每个contig
        assembled_contigs = []
        for i, contig_nodes in enumerate(contigs):
            logger.info(f"正在组装contig {i+1}/{len(contigs)}, 包含 {len(contig_nodes)} 个序列")
            
            # 组装序列
            assembled_seq = self.assemble_contig(contig_nodes)
            
            # 生成共识序列
            contig_sequences = [self.sequences[node] for node in contig_nodes]
            contig_qualities = [self.qualities[node] for node in contig_nodes]
            consensus_seq = self.generate_consensus(contig_sequences, contig_qualities)
            
            assembled_contigs.append({
                'contig_id': f"contig_{i+1}",
                'nodes': contig_nodes,
                'assembled_sequence': assembled_seq,
                'consensus_sequence': consensus_seq,
                'length': len(consensus_seq),
                'num_sequences': len(contig_nodes)
            })
            
            logger.info(f"Contig {i+1} 组装完成，长度: {len(consensus_seq)}")
        
        self.contigs = assembled_contigs
        return assembled_contigs
    
    def evaluate_assembly_quality(self, assembled_seq: str, qualities: List[int], 
                                 overlap_length: int, forward_length: int, reverse_length: int) -> Dict:
        """评估组装质量"""
        quality_metrics = {}
        
        # 1. 平均质量分数
        if qualities:
            quality_metrics['average_quality'] = np.mean(qualities)
            quality_metrics['min_quality'] = np.min(qualities)
            quality_metrics['max_quality'] = np.max(qualities)
        else:
            quality_metrics['average_quality'] = 0
            quality_metrics['min_quality'] = 0
            quality_metrics['max_quality'] = 0
        
        # 2. 序列长度评估
        quality_metrics['sequence_length'] = len(assembled_seq)
        quality_metrics['coverage_ratio'] = (forward_length + reverse_length) / len(assembled_seq) if len(assembled_seq) > 0 else 0
        
        # 3. 重叠区域质量评估
        if overlap_length > 0:
            quality_metrics['overlap_ratio'] = overlap_length / len(assembled_seq)
            quality_metrics['overlap_quality'] = np.mean(qualities[:overlap_length]) if qualities and len(qualities) >= overlap_length else 0
        else:
            quality_metrics['overlap_ratio'] = 0
            quality_metrics['overlap_quality'] = 0
        
        # 4. GC含量计算
        gc_count = assembled_seq.count('G') + assembled_seq.count('C')
        quality_metrics['gc_content'] = gc_count / len(assembled_seq) if len(assembled_seq) > 0 else 0
        
        # 5. 低质量碱基比例
        if qualities:
            low_quality_bases = sum(1 for q in qualities if q < 20)
            quality_metrics['low_quality_ratio'] = low_quality_bases / len(qualities)
        else:
            quality_metrics['low_quality_ratio'] = 0
        
        # 6. 组装质量评分 (0-100)
        quality_score = 0
        
        # 质量分数贡献 (40%)
        if quality_metrics['average_quality'] >= 30:
            quality_score += 40
        elif quality_metrics['average_quality'] >= 20:
            quality_score += 30
        elif quality_metrics['average_quality'] >= 10:
            quality_score += 20
        else:
            quality_score += 10
        
        # 重叠比例贡献 (30%)
        if quality_metrics['overlap_ratio'] >= 0.1:
            quality_score += 30
        elif quality_metrics['overlap_ratio'] >= 0.05:
            quality_score += 20
        elif quality_metrics['overlap_ratio'] >= 0.02:
            quality_score += 10
        
        # 序列长度贡献 (20%)
        if quality_metrics['sequence_length'] >= 1000:
            quality_score += 20
        elif quality_metrics['sequence_length'] >= 500:
            quality_score += 15
        elif quality_metrics['sequence_length'] >= 200:
            quality_score += 10
        
        # 低质量碱基比例贡献 (10%)
        if quality_metrics['low_quality_ratio'] <= 0.1:
            quality_score += 10
        elif quality_metrics['low_quality_ratio'] <= 0.2:
            quality_score += 5
        
        quality_metrics['overall_quality_score'] = quality_score
        
        # 7. 质量等级评估
        if quality_score >= 80:
            quality_metrics['quality_grade'] = 'A'
        elif quality_score >= 60:
            quality_metrics['quality_grade'] = 'B'
        elif quality_score >= 40:
            quality_metrics['quality_grade'] = 'C'
        else:
            quality_metrics['quality_grade'] = 'D'
        
        return quality_metrics
    
    def save_results(self, output_dir: str = "sanger_assembly_results"):
        """保存组装结果"""
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # 保存FASTA文件
        fasta_file = os.path.join(output_dir, "assembled_sequences.fasta")
        with open(fasta_file, 'w') as f:
            for sample in self.assembled_samples:
                record = SeqRecord(
                    Seq(sample['assembled_sequence']),
                    id=sample['sample_name'],
                    description=f"Assembled sequence, length={sample['length']}, overlap={sample['overlap_length']}, quality_score={sample.get('quality_score', 0)}"
                )
                SeqIO.write(record, f, "fasta")
        
        # 保存结果摘要
        summary_file = os.path.join(output_dir, "assembly_summary.csv")
        summary_data = []
        for sample in self.assembled_samples:
            summary_data.append({
                'Sample': sample['sample_name'],
                'Length': sample['length'],
                'Overlap_Length': sample['overlap_length'],
                'Forward_Length': sample['forward_length'],
                'Reverse_Length': sample['reverse_length'],
                'Average_Quality': np.mean(sample['qualities']) if sample['qualities'] else 0,
                'Quality_Score': sample.get('quality_score', 0),
                'Quality_Grade': sample.get('quality_grade', 'N/A'),
                'GC_Content': sample.get('gc_content', 0),
                'Low_Quality_Ratio': sample.get('low_quality_ratio', 0),
                'Overlap_Ratio': sample.get('overlap_ratio', 0)
            })
        
        df = pd.DataFrame(summary_data)
        df.to_csv(summary_file, index=False)
        
        # 保存质量评估报告
        quality_report_file = os.path.join(output_dir, "quality_report.txt")
        with open(quality_report_file, 'w', encoding='utf-8') as f:
            f.write("Sanger测序组装质量评估报告\n")
            f.write("=" * 50 + "\n\n")
            
            # 总体统计
            total_samples = len(self.assembled_samples)
            avg_quality_score = np.mean([s.get('quality_score', 0) for s in self.assembled_samples])
            avg_length = np.mean([s['length'] for s in self.assembled_samples])
            avg_overlap = np.mean([s['overlap_length'] for s in self.assembled_samples])
            
            f.write(f"总体统计:\n")
            f.write(f"  成功组装样本数: {total_samples}\n")
            f.write(f"  平均质量评分: {avg_quality_score:.2f}\n")
            f.write(f"  平均序列长度: {avg_length:.0f} bp\n")
            f.write(f"  平均重叠长度: {avg_overlap:.0f} bp\n\n")
            
            # 质量等级分布
            grade_counts = {}
            for sample in self.assembled_samples:
                grade = sample.get('quality_grade', 'N/A')
                grade_counts[grade] = grade_counts.get(grade, 0) + 1
            
            f.write("质量等级分布:\n")
            for grade in ['A', 'B', 'C', 'D']:
                count = grade_counts.get(grade, 0)
                percentage = (count / total_samples * 100) if total_samples > 0 else 0
                f.write(f"  {grade}级: {count}个样本 ({percentage:.1f}%)\n")
            f.write("\n")
            
            # 详细样本信息
            f.write("详细样本信息:\n")
            f.write("-" * 80 + "\n")
            for sample in self.assembled_samples:
                f.write(f"样本: {sample['sample_name']}\n")
                f.write(f"  序列长度: {sample['length']} bp\n")
                f.write(f"  重叠长度: {sample['overlap_length']} bp\n")
                f.write(f"  平均质量: {np.mean(sample['qualities']):.2f}\n")
                f.write(f"  质量评分: {sample.get('quality_score', 0)}/100\n")
                f.write(f"  质量等级: {sample.get('quality_grade', 'N/A')}\n")
                f.write(f"  GC含量: {sample.get('gc_content', 0):.3f}\n")
                f.write(f"  低质量碱基比例: {sample.get('low_quality_ratio', 0):.3f}\n")
                f.write("-" * 40 + "\n")
        
        logger.info(f"结果已保存到 {output_dir}")
        logger.info(f"质量评估报告: {quality_report_file}")
        return fasta_file, summary_file, quality_report_file

def main():
    """主函数"""
    # 初始化组装器
    assembler = GraphBasedAssembler(
        min_overlap=50,
        min_quality=20,
        min_length=100
    )
    
    # 获取当前工作目录
    data_dir = os.getcwd()
    logger.info(f"数据目录: {data_dir}")
    
    # 读取ab1文件
    samples = assembler.read_ab1_files(data_dir)
    
    if not samples:
        logger.error("未找到任何ab1文件")
        return
    
    # 处理所有样本
    assembled_samples = assembler.process_samples(samples, data_dir)
    
    if assembled_samples:
        # 保存结果
        fasta_file, summary_file, quality_report_file = assembler.save_results()
        logger.info(f"组装完成！")
        logger.info(f"FASTA文件: {fasta_file}")
        logger.info(f"摘要文件: {summary_file}")
        logger.info(f"质量评估报告: {quality_report_file}")
        logger.info(f"共成功组装 {len(assembled_samples)} 个样本")
        
        # 显示质量统计
        avg_quality_score = np.mean([s.get('quality_score', 0) for s in assembled_samples])
        logger.info(f"平均质量评分: {avg_quality_score:.2f}/100")
    else:
        logger.error("没有成功组装的样本")

if __name__ == "__main__":
    main() 