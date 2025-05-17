import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv
import sys
import os

# probe合成用データセット
initiator_seq = {
    'S23': 'GGGTGGTCGTCGAAGTCGTAT',
    'S41': 'GCTCGACGTTCCTTTGCAACA',
    'S45': 'CCTCCACGTTCCATCTAAGCT',
    'S72': 'CGGTGGAGTGGCAAGTAGGAT',
    'S73': 'CGGTCAGGTGGCTAGTATGGA',
    'A161': 'GGTACGCGAAGGTAGGTGTAA'
}


def write_output_file(probe_list, output_directory, target_name, hairpin, output_format="csv"):
    output_path = os.path.join(output_directory, f'HCR_{target_name}_{hairpin}_Probe.{output_format}')
    
    if output_format == 'fasta':
        probe_fasta = [SeqRecord(Seq(seq_f), id=name_f) for name_f, seq_f in probe_list.items()]
        SeqIO.write(probe_fasta, output_path, 'fasta')
    elif output_format == 'csv':
        with open(output_path, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['OligoName', 'Sequence'])
            for name_c, seq_c in probe_list.items():
                writer.writerow([name_c, seq_c])
    else:
        print('指定できない出力フォーマットです')
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", required=True, help="Input file path")
    parser.add_argument("-t", "--target_name", required=True, help="Target name")
    parser.add_argument("-H", "--hairpin", required=True, help="Hairpin No.")
    parser.add_argument("-f", "--output_format", choices=["fasta", "csv"], default="csv", help="Output format")
    parser.add_argument("-o", "--output_directory", default=".", help="Output directory")
    args = parser.parse_args()

    # probe_regionファイル読み込み
    try:
        records = SeqIO.parse(args.input_file, 'fasta')
    except FileNotFoundError:
        print(f"ファイルが見つかりません: {args.input_file}")
        sys.exit(1)

    # ヘアピンDNA numberの確認
    if args.hairpin not in initiator_seq:
        print('存在しないヘアピンDNA numberです')
        sys.exit(1)

    # probeを作成
    probe_list = {}
    for record in records:
        P1 = record.seq[:25]
        P2 = record.seq[27:]
        P1_probe = initiator_seq[args.hairpin][:9] + 'AA' + P1.reverse_complement().lower()
        P2_probe = P2.reverse_complement().lower() + 'AA' + initiator_seq[args.hairpin][9:]

        probe_list[f'{record.name}_{args.hairpin}_P1'] = str(P1_probe)
        probe_list[f'{record.name}_{args.hairpin}_P2'] = str(P2_probe)

    write_output_file(probe_list, args.output_directory, args.target_name, args.hairpin, args.output_format)

if __name__ == "__main__":
    main()
