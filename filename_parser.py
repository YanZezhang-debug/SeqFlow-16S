#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
智能文件名解析器
支持多种命名规则的自动识别
"""

import re
import os
from typing import Dict, List, Tuple, Optional
from collections import defaultdict
import yaml


class FilenameParser:
    """智能文件名解析器"""
    
    def __init__(self, config_file: str = "config.yaml"):
        """
        初始化解析器
        
        参数:
            config_file: 配置文件路径
        """
        self.patterns = []
        self.load_config(config_file)
    
    def load_config(self, config_file: str):
        """加载配置文件"""
        if not os.path.exists(config_file):
            print(f"警告: 配置文件 {config_file} 不存在，使用默认规则")
            self._load_default_patterns()
            return
        
        try:
            with open(config_file, 'r', encoding='utf-8') as f:
                config = yaml.safe_load(f)
            
            if 'filename_patterns' in config:
                self.patterns = config['filename_patterns']
            else:
                print("警告: 配置文件中未找到 filename_patterns，使用默认规则")
                self._load_default_patterns()
        except Exception as e:
            print(f"警告: 加载配置文件失败 ({e})，使用默认规则")
            self._load_default_patterns()
    
    def _load_default_patterns(self):
        """加载默认的文件名规则"""
        self.patterns = [
            {
                'name': '标准格式_括号样本名',
                'pattern': r'^\d+_\d+_\(([^)]+)\)_\[16S([FR])\]',
                'sample_group': 1,
                'direction_group': 2,
                'description': '格式: 编号_订单号_(样本ID)_[16S方向]'
            },
            {
                'name': '简化格式_下划线分隔',
                'pattern': r'^([^_]+)_([FR])\.ab1$',
                'sample_group': 1,
                'direction_group': 2,
                'description': '格式: 样本ID_方向'
            },
            {
                'name': '简化格式_横杠分隔',
                'pattern': r'^([^-]+)-([FR])\.ab1$',
                'sample_group': 1,
                'direction_group': 2,
                'description': '格式: 样本ID-方向'
            }
        ]
    
    def parse_filename(self, filename: str) -> Optional[Dict[str, str]]:
        """
        解析单个文件名
        
        参数:
            filename: 文件名
        
        返回:
            包含样本ID、方向、匹配规则等信息的字典，解析失败返回None
        """
        basename = os.path.basename(filename)
        
        # 尝试所有规则
        for pattern_info in self.patterns:
            pattern = pattern_info['pattern']
            match = re.search(pattern, basename)
            
            if match:
                sample_id = match.group(pattern_info['sample_group'])
                direction = match.group(pattern_info['direction_group'])
                
                return {
                    'filename': filename,
                    'basename': basename,
                    'sample_id': sample_id,
                    'direction': direction,
                    'pattern_name': pattern_info['name'],
                    'pattern_description': pattern_info['description']
                }
        
        return None
    
    def parse_directory(self, directory: str, extension: str = ".ab1") -> Dict[str, List[Dict]]:
        """
        解析目录中的所有文件
        
        参数:
            directory: 目录路径
            extension: 文件扩展名
        
        返回:
            按样本ID分组的文件信息字典
        """
        # 获取所有指定扩展名的文件
        files = [f for f in os.listdir(directory) if f.endswith(extension)]
        
        # 解析每个文件
        parsed_files = []
        unparsed_files = []
        
        for filename in files:
            filepath = os.path.join(directory, filename)
            result = self.parse_filename(filename)
            
            if result:
                result['filepath'] = filepath
                parsed_files.append(result)
            else:
                unparsed_files.append(filename)
        
        # 按样本ID分组
        grouped = defaultdict(list)
        for file_info in parsed_files:
            grouped[file_info['sample_id']].append(file_info)
        
        return {
            'grouped': dict(grouped),
            'parsed_files': parsed_files,
            'unparsed_files': unparsed_files,
            'total_files': len(files),
            'parsed_count': len(parsed_files),
            'unparsed_count': len(unparsed_files)
        }
    
    def find_pairs(self, directory: str, extension: str = ".ab1") -> List[Dict]:
        """
        查找配对的正反向测序文件
        
        参数:
            directory: 目录路径
            extension: 文件扩展名
        
        返回:
            配对信息列表
        """
        result = self.parse_directory(directory, extension)
        grouped = result['grouped']
        
        pairs = []
        singles = []
        
        for sample_id, files in grouped.items():
            # 按方向分组
            by_direction = {'F': [], 'R': []}
            for file_info in files:
                by_direction[file_info['direction']].append(file_info)
            
            # 检查是否有配对
            if by_direction['F'] and by_direction['R']:
                # 如果有多个，取第一个（可以根据需要调整）
                pairs.append({
                    'sample_id': sample_id,
                    'forward': by_direction['F'][0],
                    'reverse': by_direction['R'][0],
                    'has_pair': True
                })
            else:
                # 单个文件
                if by_direction['F']:
                    singles.append({
                        'sample_id': sample_id,
                        'file': by_direction['F'][0],
                        'direction': 'F',
                        'has_pair': False
                    })
                if by_direction['R']:
                    singles.append({
                        'sample_id': sample_id,
                        'file': by_direction['R'][0],
                        'direction': 'R',
                        'has_pair': False
                    })
        
        return {
            'pairs': pairs,
            'singles': singles,
            'unparsed_files': result['unparsed_files'],
            'statistics': {
                'total_files': result['total_files'],
                'parsed_files': result['parsed_count'],
                'unparsed_files': result['unparsed_count'],
                'paired_samples': len(pairs),
                'single_files': len(singles)
            }
        }
    
    def print_summary(self, directory: str, extension: str = ".ab1"):
        """
        打印解析摘要
        
        参数:
            directory: 目录路径
            extension: 文件扩展名
        """
        result = self.find_pairs(directory, extension)
        stats = result['statistics']
        
        print("\n" + "=" * 80)
        print("文件名解析摘要")
        print("=" * 80)
        print(f"目录: {directory}")
        print(f"文件扩展名: {extension}")
        print(f"\n总文件数: {stats['total_files']}")
        print(f"  ├─ 成功解析: {stats['parsed_files']}")
        print(f"  └─ 无法解析: {stats['unparsed_files']}")
        print(f"\n样本统计:")
        print(f"  ├─ 配对样本: {stats['paired_samples']}")
        print(f"  └─ 单个文件: {stats['single_files']}")
        
        if result['pairs']:
            print(f"\n配对样本列表:")
            for i, pair in enumerate(result['pairs'], 1):
                print(f"  {i}. {pair['sample_id']}")
                print(f"     F: {pair['forward']['basename']}")
                print(f"     R: {pair['reverse']['basename']}")
        
        if result['singles']:
            print(f"\n单个文件列表:")
            for i, single in enumerate(result['singles'], 1):
                print(f"  {i}. {single['sample_id']} ({single['direction']}): {single['file']['basename']}")
        
        if result['unparsed_files']:
            print(f"\n无法解析的文件:")
            for i, filename in enumerate(result['unparsed_files'], 1):
                print(f"  {i}. {filename}")
            print(f"\n提示: 这些文件的命名格式不符合任何已知规则")
            print(f"      可以在 config.yaml 中添加自定义规则")
        
        print("=" * 80 + "\n")
        
        return result
    
    def get_supported_patterns(self) -> List[Dict]:
        """获取支持的所有命名规则"""
        return self.patterns
    
    def print_supported_patterns(self):
        """打印支持的命名规则"""
        print("\n" + "=" * 80)
        print("支持的文件命名规则")
        print("=" * 80)
        
        for i, pattern in enumerate(self.patterns, 1):
            print(f"\n{i}. {pattern['name']}")
            print(f"   描述: {pattern['description']}")
            print(f"   正则: {pattern['pattern']}")
        
        print("\n" + "=" * 80)
        print("如需添加自定义规则，请编辑 config.yaml 文件")
        print("=" * 80 + "\n")


def main():
    """测试函数"""
    import sys
    
    parser = FilenameParser()
    
    if len(sys.argv) > 1:
        if sys.argv[1] == '--patterns':
            parser.print_supported_patterns()
        else:
            directory = sys.argv[1]
            if os.path.isdir(directory):
                parser.print_summary(directory)
            else:
                print(f"错误: {directory} 不是有效的目录")
    else:
        # 使用当前目录
        parser.print_summary(".")


if __name__ == "__main__":
    main()

