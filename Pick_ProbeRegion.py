import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction
import pandas as pd
import sys
import os

# GC含有率の条件を満たし、TTのない領域のみ取り出す関数
def pick_probe_region(record, target_name, GC_min=45, GC_max=55):
    region_list = []
    
    for i in range(len(record.seq) - 51):
        sub_sequence = record.seq[i:i + 52].upper()
        P1 = sub_sequence[:25]
        P2 = sub_sequence[27:]
        P1_GC = gc_fraction(P1) * 100
        P2_GC = gc_fraction(P2) * 100
        inter_bindsite = sub_sequence[25:27]

        if (GC_min <= P1_GC <= GC_max) and (GC_min <= P2_GC <= GC_max) and (inter_bindsite != 'TT'):
            region_list.append([
                f'{target_name}_Region{i + 1}', sub_sequence, i + 1, i + 52, 
                P1, P2, P1_GC, P2_GC
            ])
    
    return pd.DataFrame(region_list, columns=[
        'region_id', 'sequence', 'start base', 'end base', 'P1', 'P2', 'P1 GC fraction', 'P2 GC fraction'
    ]).set_index('region_id')

# 基準となる領域を決めて重複している領域を除く関数
def remove_overlap(region_df, stem_region_num):
    # あらかじめソートしておく
    region_df = region_df.sort_values('start base').copy()
    region_df = region_df.reset_index()
    
    stem_region = region_df.iloc[stem_region_num]
    new_index = [stem_region_num]

    # 後方チェック
    for i in range(stem_region_num + 1, len(region_df)):
        if region_df.loc[i, 'start base'] > stem_region['end base']:
            new_index.append(i)
            stem_region = region_df.loc[i]

    # 前方チェック（逆順）
    stem_region = region_df.iloc[stem_region_num]
    for i in range(stem_region_num - 1, -1, -1):
        if region_df.loc[i, 'end base'] < stem_region['start base']:
            new_index.append(i)
            stem_region = region_df.loc[i]

    filtered_df = region_df.loc[new_index]
    return filtered_df.sort_values('start base').set_index('region_id')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", required=True, help="Input file path")
    parser.add_argument("-m", "--gc_min", default=45, help="GC min")
    parser.add_argument("-M", "--gc_max", default=55, help="GC max")
    parser.add_argument("-t", "--target_name", required=True, help="target name")
    parser.add_argument("-o", "--output_directory", required=True, help="Output directory")
    args = parser.parse_args()

    try:
        record = SeqIO.read(args.input_file, 'fasta')
    except FileNotFoundError:
        print(f"ファイルが見つかりません: {args.input_file}")
        sys.exit(1)

    # ProbeRegionの選定
    region_df = pick_probe_region(record, args.target_name, args.gc_min, args.gc_max)

    # 最も多くの非重複領域を得る stem_region_num を探索
    max_count = -1
    best_df = None
    for i in range(len(region_df)):
        try:
            filtered = remove_overlap(region_df, i)
            if len(filtered) > max_count:
                max_count = len(filtered)
                best_df = filtered
        except Exception:
            continue  # 念のため安全にスキップ

    region_df_filtered = best_df

    # CSVファイル出力
    output_csv_path = os.path.join(args.output_directory, f'{args.target_name}_ProbeRegion.csv')
    region_df_filtered.to_csv(output_csv_path)
    
    # Blast用のfastaファイルを保存
    fasta_records = [
        SeqRecord(Seq(row['sequence']), id=row.name) for _, row in region_df_filtered.iterrows()
    ]
    output_fasta_path = os.path.join(args.output_directory, f'{args.target_name}_ProbeRegion.fasta')
    SeqIO.write(fasta_records, output_fasta_path, 'fasta')

if __name__ == "__main__":
    main()