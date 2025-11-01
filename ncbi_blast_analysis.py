#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
NCBI BLAST分析脚本
用于对组装后的16S rRNA序列进行BLAST分析并保存结果到Excel
"""

import os
import sys
import time
import logging
import json
import argparse
from typing import List, Dict, Optional
from datetime import datetime
import pandas as pd
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
import warnings
warnings.filterwarnings('ignore')

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('blast_analysis.log', encoding='utf-8'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)


class NCBIBlastAnalyzer:
    """NCBI BLAST分析器"""
    
    def __init__(self, database: str = "nt", max_target_seqs: int = 10, program: str = "blastn"):
        """
        初始化BLAST分析器
        
        参数:
            database: BLAST数据库名称
                     - "nt": 核苷酸数据库（默认，最全面但较慢）
                     - "refseq_rna": RefSeq RNA数据库（推荐用于16S）
                     - "16S_ribosomal_RNA": 16S专用数据库（如果可用）
            max_target_seqs: 返回的最大匹配数，默认10条
            program: BLAST程序类型，默认"blastn"
        """
        self.database = database
        self.max_target_seqs = max_target_seqs
        self.program = program
        self.blast_results = []
        
    def read_fasta(self, fasta_file: str) -> List[Dict]:
        """
        读取FASTA文件
        
        参数:
            fasta_file: FASTA文件路径
            
        返回:
            序列列表
        """
        logger.info(f"正在读取FASTA文件: {fasta_file}")
        sequences = []
        
        try:
            for record in SeqIO.parse(fasta_file, "fasta"):
                sequences.append({
                    'id': record.id,
                    'description': record.description,
                    'sequence': str(record.seq),
                    'length': len(record.seq)
                })
            logger.info(f"成功读取 {len(sequences)} 条序列")
            return sequences
        except Exception as e:
            logger.error(f"读取FASTA文件失败: {e}")
            return []
    
    def run_blast(self, sequence: str, sequence_id: str, retry_count: int = 3) -> Optional[str]:
        """
        对单个序列运行BLAST，支持重试机制
        
        参数:
            sequence: DNA序列
            sequence_id: 序列ID
            retry_count: 失败后的重试次数
            
        返回:
            BLAST结果XML字符串
        """
        logger.info(f"正在对序列 {sequence_id} 运行BLAST (长度: {len(sequence)} bp)...")
        
        for attempt in range(retry_count):
            try:
                # 调用NCBI BLAST API
                # 使用blastn程序，针对16S核糖体RNA数据库
                result_handle = NCBIWWW.qblast(
                    program=self.program,
                    database=self.database,
                    sequence=sequence,
                    hitlist_size=self.max_target_seqs,
                    expect=10,  # E值阈值
                    megablast=True  # 使用megablast获得更快的结果
                )
                
                # 读取结果
                blast_results = result_handle.read()
                result_handle.close()
                
                logger.info(f"序列 {sequence_id} BLAST完成")
                return blast_results
                
            except Exception as e:
                if attempt < retry_count - 1:
                    logger.warning(f"序列 {sequence_id} BLAST失败 (尝试 {attempt + 1}/{retry_count}): {e}")
                    logger.info(f"等待10秒后重试...")
                    time.sleep(10)
                else:
                    logger.error(f"序列 {sequence_id} BLAST失败，已达到最大重试次数: {e}")
                    return None
        
        return None
    
    def parse_blast_results(self, blast_xml: str, sequence_id: str) -> List[Dict]:
        """
        解析BLAST XML结果，只返回评分最高的一个匹配
        
        参数:
            blast_xml: BLAST结果XML字符串
            sequence_id: 序列ID
            
        返回:
            仅包含最佳匹配的结果列表（单个元素）
        """
        from io import BytesIO
        
        results = []
        
        try:
            # 将字符串转换为字节流
            if isinstance(blast_xml, str):
                blast_xml = blast_xml.encode('utf-8')
            
            # 解析XML - 使用read而不是parse（因为只有一个查询）
            blast_record = NCBIXML.read(BytesIO(blast_xml))
            
            # 只取第一个比对结果（BLAST默认按评分从高到低排序）
            if blast_record.alignments:
                alignment = blast_record.alignments[0]  # 评分最高的匹配
                hsp = alignment.hsps[0]  # 最佳的高分段配对
                
                # 计算相似度百分比
                identity_percent = (hsp.identities / hsp.align_length) * 100
                coverage_percent = (hsp.align_length / blast_record.query_length) * 100
                
                # 提取物种信息
                hit_def = alignment.hit_def
                # 尝试提取菌名（通常在描述的开头）
                species_name = hit_def.split(',')[0] if ',' in hit_def else hit_def
                
                result = {
                    '样本ID': sequence_id,
                    '查询长度': blast_record.query_length,
                    '匹配序列ID': alignment.accession,
                    '物种名称': species_name,
                    '完整描述': hit_def,
                    '比对长度': hsp.align_length,
                    '相似碱基数': hsp.identities,
                    '相似度(%)': round(identity_percent, 2),
                    '覆盖度(%)': round(coverage_percent, 2),
                    'E值': hsp.expect,
                    '比对得分': hsp.score,
                    '比特分数': hsp.bits,
                    '查询起始位置': hsp.query_start,
                    '查询结束位置': hsp.query_end,
                    '匹配起始位置': hsp.sbjct_start,
                    '匹配结束位置': hsp.sbjct_end,
                    '正链/负链': '+' if hsp.sbjct_start < hsp.sbjct_end else '-'
                }
                
                results.append(result)
                logger.info(f"序列 {sequence_id} 最佳匹配: {species_name} (相似度: {identity_percent:.2f}%)")
            else:
                logger.warning(f"序列 {sequence_id} 没有找到BLAST匹配")
                
        except Exception as e:
            logger.error(f"解析BLAST结果失败: {e}")
            import traceback
            logger.error(traceback.format_exc())
        
        return results
    
    def save_progress(self, progress_file: str = "blast_progress.json"):
        """保存当前进度"""
        try:
            with open(progress_file, 'w', encoding='utf-8') as f:
                json.dump(self.blast_results, f, ensure_ascii=False, indent=2)
            logger.info(f"进度已保存到 {progress_file}")
        except Exception as e:
            logger.error(f"保存进度失败: {e}")
    
    def load_progress(self, progress_file: str = "blast_progress.json") -> List[str]:
        """
        加载之前的进度
        
        返回:
            已完成的序列ID列表
        """
        if not os.path.exists(progress_file):
            return []
        
        try:
            with open(progress_file, 'r', encoding='utf-8') as f:
                self.blast_results = json.load(f)
            completed_ids = [r['样本ID'] for r in self.blast_results]
            logger.info(f"从 {progress_file} 加载了 {len(completed_ids)} 条已完成的结果")
            return completed_ids
        except Exception as e:
            logger.error(f"加载进度失败: {e}")
            return []
    
    def analyze_sequences(self, fasta_file: str, delay: int = 10, resume: bool = True):
        """
        分析FASTA文件中的所有序列
        
        参数:
            fasta_file: FASTA文件路径
            delay: 每次BLAST之间的延迟（秒），避免请求过于频繁
            resume: 是否从上次中断的地方继续
        """
        # 读取序列
        sequences = self.read_fasta(fasta_file)
        
        if not sequences:
            logger.error("未找到任何序列")
            return
        
        # 加载之前的进度
        completed_ids = []
        if resume:
            completed_ids = self.load_progress()
        
        # 过滤掉已完成的序列
        remaining_sequences = [s for s in sequences if s['id'] not in completed_ids]
        
        if not remaining_sequences:
            logger.info("所有序列都已完成BLAST分析！")
            return
        
        logger.info(f"开始分析 {len(remaining_sequences)} 条序列 (总共 {len(sequences)} 条，已完成 {len(completed_ids)} 条)...")
        
        # 对每条序列运行BLAST
        for idx, seq_data in enumerate(remaining_sequences, 1):
            logger.info(f"\n{'='*60}")
            logger.info(f"进度: {idx}/{len(remaining_sequences)} (总进度: {len(completed_ids) + idx}/{len(sequences)})")
            logger.info(f"正在处理序列: {seq_data['id']}")
            
            # 运行BLAST
            blast_xml = self.run_blast(seq_data['sequence'], seq_data['id'])
            
            if blast_xml:
                # 解析结果
                results = self.parse_blast_results(blast_xml, seq_data['id'])
                self.blast_results.extend(results)
                
                # 保存进度
                self.save_progress()
            
            # 如果不是最后一条序列，等待一段时间再继续
            if idx < len(remaining_sequences):
                logger.info(f"等待 {delay} 秒后继续下一个序列...")
                time.sleep(delay)
        
        logger.info(f"\n所有序列分析完成！共获得 {len(self.blast_results)} 条BLAST结果")
    
    def save_to_excel(self, output_file: str = "blast_results.xlsx"):
        """
        保存结果到Excel文件
        
        参数:
            output_file: 输出Excel文件路径
        """
        if not self.blast_results:
            logger.warning("没有BLAST结果可保存")
            return
        
        logger.info(f"正在保存结果到 {output_file}...")
        
        # 创建DataFrame（每个序列只有一个最佳匹配）
        df = pd.DataFrame(self.blast_results)
        
        # 按相似度排序
        df_sorted = df.sort_values('相似度(%)', ascending=False)
        
        # 保存到Excel
        df_sorted.to_excel(output_file, index=False, engine='openpyxl')
        
        logger.info(f"结果已成功保存到 {output_file}")
        logger.info(f"  - 共 {len(df_sorted)} 个序列的最佳BLAST匹配结果")
    
    def generate_report(self, output_file: str = "blast_analysis_report.txt"):
        """
        生成文本格式的分析报告
        
        参数:
            output_file: 输出报告文件路径
        """
        if not self.blast_results:
            logger.warning("没有BLAST结果可生成报告")
            return
        
        logger.info(f"正在生成分析报告: {output_file}")
        
        df = pd.DataFrame(self.blast_results)
        # 按相似度排序
        df_sorted = df.sort_values('相似度(%)', ascending=False)
        
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write("NCBI BLAST分析报告\n")
            f.write("=" * 80 + "\n\n")
            f.write(f"分析时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"数据库: {self.database}\n")
            f.write(f"分析序列数: {len(df_sorted)}\n")
            f.write(f"(每个序列只显示评分最高的匹配)\n\n")
            
            f.write("=" * 80 + "\n")
            f.write("各序列最佳匹配结果\n")
            f.write("=" * 80 + "\n\n")
            
            for idx, row in df_sorted.iterrows():
                f.write(f"序列: {row['样本ID']}\n")
                f.write(f"  查询序列长度: {row['查询长度']} bp\n")
                f.write(f"  最佳匹配物种: {row['物种名称']}\n")
                f.write(f"  GenBank登录号: {row['匹配序列ID']}\n")
                f.write(f"  相似度: {row['相似度(%)']}%\n")
                f.write(f"  覆盖度: {row['覆盖度(%)']}%\n")
                f.write(f"  E值: {row['E值']}\n")
                f.write(f"  比特分数: {row['比特分数']}\n")
                f.write("-" * 80 + "\n")
        
        logger.info(f"报告已保存到 {output_file}")


def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='NCBI BLAST分析脚本 - 对16S rRNA序列进行在线BLAST分析',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 使用默认设置运行
  python ncbi_blast_analysis.py
  
  # 指定输入文件
  python ncbi_blast_analysis.py -i my_sequences.fasta
  
  # 使用RefSeq数据库（更快）
  python ncbi_blast_analysis.py -d refseq_rna
  
  # 调整延迟时间
  python ncbi_blast_analysis.py --delay 15
  
  # 不使用断点续传（从头开始）
  python ncbi_blast_analysis.py --no-resume
        """
    )
    
    parser.add_argument(
        '-i', '--input',
        default='sanger_assembly_results/assembled_sequences.fasta',
        help='输入FASTA文件路径 (默认: sanger_assembly_results/assembled_sequences.fasta)'
    )
    
    parser.add_argument(
        '-d', '--database',
        default='nt',
        choices=['nt', 'refseq_rna', '16S_ribosomal_RNA'],
        help='BLAST数据库选择 (默认: nt)'
    )
    
    parser.add_argument(
        '-n', '--max-hits',
        type=int,
        default=10,
        help='每个序列返回的最大匹配数 (默认: 10)'
    )
    
    parser.add_argument(
        '--delay',
        type=int,
        default=10,
        help='每次BLAST之间的延迟秒数 (默认: 10)'
    )
    
    parser.add_argument(
        '--no-resume',
        action='store_true',
        help='不使用断点续传，从头开始分析'
    )
    
    parser.add_argument(
        '-o', '--output-dir',
        default='blast_results',
        help='输出目录 (默认: blast_results)'
    )
    
    return parser.parse_args()


def main():
    """主函数"""
    # 解析命令行参数
    args = parse_arguments()
    
    # 检查输入文件是否存在
    if not os.path.exists(args.input):
        logger.error(f"未找到FASTA文件: {args.input}")
        logger.error("请先运行 graph_based_assembly.py 进行序列组装")
        logger.error("或使用 -i 参数指定正确的FASTA文件路径")
        return
    
    # 创建BLAST分析器
    logger.info("=" * 80)
    logger.info("NCBI BLAST分析器配置")
    logger.info("=" * 80)
    logger.info(f"输入文件: {args.input}")
    logger.info(f"数据库: {args.database}")
    logger.info(f"最大匹配数: {args.max_hits}")
    logger.info(f"延迟时间: {args.delay} 秒")
    logger.info(f"断点续传: {'否' if args.no_resume else '是'}")
    logger.info(f"输出目录: {args.output_dir}")
    logger.info("=" * 80 + "\n")
    
    analyzer = NCBIBlastAnalyzer(
        database=args.database,
        max_target_seqs=args.max_hits
    )
    
    # 分析序列
    logger.info("\n" + "=" * 80)
    logger.info("开始BLAST分析")
    logger.info("=" * 80 + "\n")
    logger.info("注意: NCBI BLAST是在线服务，每个序列大约需要20-60秒")
    logger.info("如果中断，可以重新运行脚本继续分析（断点续传）")
    logger.info("请耐心等待...\n")
    
    try:
        analyzer.analyze_sequences(
            fasta_file=args.input,
            delay=args.delay,
            resume=not args.no_resume
        )
    except KeyboardInterrupt:
        logger.warning("\n用户中断了程序")
        logger.info("进度已保存，可以重新运行脚本继续分析")
        return
    
    # 保存结果
    if analyzer.blast_results:
        # 创建输出目录
        if not os.path.exists(args.output_dir):
            os.makedirs(args.output_dir)
        
        # 保存到Excel
        excel_file = os.path.join(args.output_dir, f"blast_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.xlsx")
        analyzer.save_to_excel(excel_file)
        
        # 生成报告
        report_file = os.path.join(args.output_dir, f"blast_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt")
        analyzer.generate_report(report_file)
        
        logger.info("\n" + "=" * 80)
        logger.info("BLAST分析完成！")
        logger.info("=" * 80)
        logger.info(f"Excel结果文件: {excel_file}")
        logger.info(f"文本报告文件: {report_file}")
        logger.info(f"日志文件: blast_analysis.log")
    else:
        logger.error("没有获得任何BLAST结果")


if __name__ == "__main__":
    main()

