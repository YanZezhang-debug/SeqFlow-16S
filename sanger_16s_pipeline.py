#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sanger 16S rRNA 测序分析完整流程
一键完成从测序文件到BLAST分析的全流程
"""

import os
import sys
import argparse
import yaml
from datetime import datetime
from typing import Dict, List
import logging

# 导入自定义模块
from filename_parser import FilenameParser


class SangerPipeline:
    """Sanger测序分析流程管理器"""
    
    def __init__(self, config_file: str = "config.yaml"):
        """
        初始化流程管理器
        
        参数:
            config_file: 配置文件路径
        """
        self.config_file = config_file
        self.config = self.load_config()
        self.setup_logging()
        self.parser = FilenameParser(config_file)
    
    def load_config(self) -> Dict:
        """加载配置文件"""
        if not os.path.exists(self.config_file):
            print(f"警告: 配置文件 {self.config_file} 不存在，使用默认配置")
            return self.get_default_config()
        
        try:
            with open(self.config_file, 'r', encoding='utf-8') as f:
                config = yaml.safe_load(f)
            return config
        except Exception as e:
            print(f"警告: 加载配置文件失败 ({e})，使用默认配置")
            return self.get_default_config()
    
    def get_default_config(self) -> Dict:
        """获取默认配置"""
        return {
            'assembly': {
                'min_overlap': 20,
                'min_quality': 20,
                'trim_ends': True
            },
            'blast': {
                'database': 'nt',
                'max_hits': 10,
                'delay': 10
            },
            'output': {
                'assembly_dir': 'sanger_assembly_results',
                'blast_dir': 'blast_results'
            },
            'logging': {
                'level': 'INFO',
                'save_to_file': True,
                'log_file': 'analysis.log'
            }
        }
    
    def setup_logging(self):
        """设置日志"""
        log_config = self.config.get('logging', {})
        level = getattr(logging, log_config.get('level', 'INFO'))
        
        handlers = [logging.StreamHandler()]
        if log_config.get('save_to_file', True):
            log_file = log_config.get('log_file', 'analysis.log')
            handlers.append(logging.FileHandler(log_file, encoding='utf-8'))
        
        logging.basicConfig(
            level=level,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=handlers
        )
        
        self.logger = logging.getLogger(__name__)
    
    def interactive_mode(self):
        """交互式模式"""
        print("\n" + "=" * 80)
        print("Sanger 16S rRNA 测序分析工具")
        print("=" * 80)
        print("\n欢迎使用！本工具将帮助您完成从测序文件到BLAST分析的全流程\n")
        
        # 1. 选择输入目录
        print("步骤 1: 选择输入目录")
        print("-" * 80)
        input_dir = input("请输入测序文件所在目录 (直接回车使用当前目录): ").strip()
        if not input_dir:
            input_dir = "."
        
        if not os.path.isdir(input_dir):
            print(f"错误: {input_dir} 不是有效的目录")
            return
        
        # 2. 扫描和解析文件
        print(f"\n步骤 2: 扫描文件")
        print("-" * 80)
        result = self.parser.print_summary(input_dir)
        
        if result['statistics']['parsed_files'] == 0:
            print("\n错误: 没有找到可识别的测序文件")
            print("提示: 请检查文件命名格式，或在 config.yaml 中添加自定义规则")
            
            view_patterns = input("\n是否查看支持的命名规则？(y/n): ").strip().lower()
            if view_patterns == 'y':
                self.parser.print_supported_patterns()
            return
        
        # 3. 确认继续
        print(f"\n步骤 3: 确认分析")
        print("-" * 80)
        print(f"找到 {result['statistics']['paired_samples']} 个配对样本")
        print(f"找到 {result['statistics']['single_files']} 个单独文件")
        
        if result['unparsed_files']:
            print(f"\n注意: 有 {len(result['unparsed_files'])} 个文件无法识别，将被跳过")
        
        confirm = input("\n是否继续进行序列组装？(y/n): ").strip().lower()
        if confirm != 'y':
            print("已取消")
            return
        
        # 4. 序列组装
        print(f"\n步骤 4: 序列组装")
        print("-" * 80)
        self.run_assembly(input_dir)
        
        # 5. BLAST分析
        print(f"\n步骤 5: BLAST分析")
        print("-" * 80)
        
        assembly_dir = self.config['output']['assembly_dir']
        fasta_file = os.path.join(assembly_dir, "assembled_sequences.fasta")
        
        if not os.path.exists(fasta_file):
            print("错误: 未找到组装后的序列文件")
            return
        
        run_blast = input("\n是否进行NCBI BLAST分析？(y/n): ").strip().lower()
        if run_blast == 'y':
            print("\n注意: BLAST分析需要连接NCBI服务器，可能需要较长时间")
            print("      每个序列大约需要20-60秒")
            
            confirm_blast = input("确认继续？(y/n): ").strip().lower()
            if confirm_blast == 'y':
                self.run_blast(fasta_file)
        
        print("\n" + "=" * 80)
        print("分析完成！")
        print("=" * 80)
    
    def run_assembly(self, input_dir: str):
        """运行序列组装"""
        self.logger.info("开始序列组装...")
        
        # 导入组装模块
        try:
            import graph_based_assembly
            
            # 创建组装器实例
            assembler = graph_based_assembly.GraphBasedAssembler(
                min_overlap=self.config['assembly'].get('min_overlap', 20),
                min_quality=self.config['assembly'].get('min_quality', 20)
            )
            
            # 运行组装
            assembler.process_directory(input_dir)
            
            self.logger.info("序列组装完成")
            
        except Exception as e:
            self.logger.error(f"序列组装失败: {e}")
            import traceback
            traceback.print_exc()
    
    def run_blast(self, fasta_file: str):
        """运行BLAST分析"""
        self.logger.info("开始BLAST分析...")
        
        try:
            import ncbi_blast_analysis
            
            # 创建BLAST分析器
            analyzer = ncbi_blast_analysis.NCBIBlastAnalyzer(
                database=self.config['blast'].get('database', 'nt'),
                max_target_seqs=self.config['blast'].get('max_hits', 10)
            )
            
            # 运行分析
            analyzer.analyze_sequences(
                fasta_file=fasta_file,
                delay=self.config['blast'].get('delay', 10),
                resume=True
            )
            
            # 保存结果
            if analyzer.blast_results:
                output_dir = self.config['output']['blast_dir']
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                
                excel_file = os.path.join(
                    output_dir,
                    f"blast_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.xlsx"
                )
                analyzer.save_to_excel(excel_file)
                
                report_file = os.path.join(
                    output_dir,
                    f"blast_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
                )
                analyzer.generate_report(report_file)
                
                self.logger.info(f"BLAST结果已保存到: {excel_file}")
            
        except KeyboardInterrupt:
            self.logger.warning("用户中断了BLAST分析")
            self.logger.info("进度已保存，可以重新运行继续分析")
        except Exception as e:
            self.logger.error(f"BLAST分析失败: {e}")
            import traceback
            traceback.print_exc()
    
    def batch_mode(self, input_dir: str, skip_blast: bool = False):
        """批处理模式（非交互）"""
        self.logger.info("=" * 80)
        self.logger.info("开始批处理模式")
        self.logger.info("=" * 80)
        
        # 扫描文件
        result = self.parser.find_pairs(input_dir)
        self.logger.info(f"找到 {result['statistics']['paired_samples']} 个配对样本")
        
        if result['statistics']['parsed_files'] == 0:
            self.logger.error("没有找到可识别的测序文件")
            return
        
        # 序列组装
        self.run_assembly(input_dir)
        
        # BLAST分析
        if not skip_blast:
            assembly_dir = self.config['output']['assembly_dir']
            fasta_file = os.path.join(assembly_dir, "assembled_sequences.fasta")
            
            if os.path.exists(fasta_file):
                self.run_blast(fasta_file)
        
        self.logger.info("=" * 80)
        self.logger.info("批处理完成")
        self.logger.info("=" * 80)


def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='Sanger 16S rRNA 测序分析完整流程',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 交互式模式（推荐新手使用）
  python sanger_16s_pipeline.py
  
  # 批处理模式
  python sanger_16s_pipeline.py -i ./data --batch
  
  # 只进行序列组装，跳过BLAST
  python sanger_16s_pipeline.py -i ./data --batch --skip-blast
  
  # 查看支持的文件命名规则
  python sanger_16s_pipeline.py --show-patterns
  
  # 只扫描文件，不进行分析
  python sanger_16s_pipeline.py -i ./data --scan-only
        """
    )
    
    parser.add_argument(
        '-i', '--input',
        help='输入目录路径'
    )
    
    parser.add_argument(
        '-c', '--config',
        default='config.yaml',
        help='配置文件路径 (默认: config.yaml)'
    )
    
    parser.add_argument(
        '--batch',
        action='store_true',
        help='批处理模式（非交互）'
    )
    
    parser.add_argument(
        '--skip-blast',
        action='store_true',
        help='跳过BLAST分析'
    )
    
    parser.add_argument(
        '--scan-only',
        action='store_true',
        help='只扫描文件，不进行分析'
    )
    
    parser.add_argument(
        '--show-patterns',
        action='store_true',
        help='显示支持的文件命名规则'
    )
    
    return parser.parse_args()


def main():
    """主函数"""
    args = parse_arguments()
    
    # 创建流程管理器
    pipeline = SangerPipeline(args.config)
    
    # 显示命名规则
    if args.show_patterns:
        pipeline.parser.print_supported_patterns()
        return
    
    # 只扫描文件
    if args.scan_only:
        input_dir = args.input or "."
        pipeline.parser.print_summary(input_dir)
        return
    
    # 批处理模式
    if args.batch:
        if not args.input:
            print("错误: 批处理模式需要指定输入目录 (-i)")
            return
        pipeline.batch_mode(args.input, args.skip_blast)
        return
    
    # 交互式模式
    pipeline.interactive_mode()


if __name__ == "__main__":
    main()

