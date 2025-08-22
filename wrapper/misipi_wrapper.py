#!/usr/bin/env python

import argparse
from pathlib import Path
import subprocess
import sys

## Extract library size from BAM file with RSamtools
def get_bam_alignment_count(bam_path):
    r_code = f"""
    library(Rsamtools)
    count <- countBam("{bam_path}")$records
    cat(count)
    """
    result = subprocess.run(["Rscript", "-e", r_code], capture_output=True, text=True, check=True)
    return int(result.stdout.strip())

## get ML table name:
def find_ml_table(directory):
    ml_files = list(directory.glob("*_ml.txt"))
    if not ml_files:
        raise FileNotFoundError(f"No *_ml.txt file found in {directory}")
    return ml_files[0].name  # Return filename only


## Main MiSiPi function
def run_misipi_main_rna(args):
    r_script = f"""
    library(MiSiPi.RNA);
    vars <- set_vars(roi = '{args.roi}',
                     bam_file = '{args.bam_file}',
                     genome = '{args.genome}',
                     plot_output = {str(args.plot_output).upper()},
                     path_to_RNAfold = '{args.path_to_RNAfold}',
                     path_to_RNAplot = '{args.path_to_RNAplot}',
                     pi_pal = '{args.pi_pal}',
                     si_pal = '{args.si_pal}',
                     annotate_region = {str(args.annotate_region).upper()},
                     weight_reads = '{args.weight_reads}',
                     gtf_file = '{args.gtf_file}',
                     write_fastas = {str(args.write_fastas).upper()},
                     out_type = '{args.out_type}');
    """

    if args.run_all:
        r_script += f"misipi_rna(vars, method = \"all\", outdir_name = \"{args.outdir_name}\");"
    else:
        if args.miRNA:
            r_script += f"misipi_rna(vars, method = \"miRNA\", outdir_name = \"{args.outdir_name}\");"
        if args.piRNA:
            r_script += f"misipi_rna(vars, method = \"piRNA\", outdir_name = \"{args.outdir_name}\");"
        if args.siRNA:
            r_script += f"misipi_rna(vars, method = \"siRNA\", outdir_name = \"{args.outdir_name}\");"

    print("Running R script with the following command:")
    print(r_script)

    # Execute the R script
    subprocess.run(["Rscript", "-e", r_script], check=True)

## Run the optional ML or HTML step
def run_postprocessing_ml_html(args):

    outdir_path = Path("misipi_output") / args.outdir_name
    if not outdir_path.exists():
        sys.exit(f"Error: Output directory {outdir_path} not found.")

    r_script = "library(MiSiPi.RNA);\n"

    if args.ml:
        print("Calculating library size from BAM...")
        library_size = get_bam_alignment_count(args.bam_file)

        # Find the full ML table path
        ml_files = list(outdir_path.glob("*_ml.txt"))
        if not ml_files:
            sys.exit(f"No *_ml.txt file found in {outdir_path}")
        ml_file = ml_files[0]
        ml_table = ml_file.name
        ml_dir = ml_file.parent.resolve()

        print(f"Using ML table: {ml_table}")
        print(f"From directory: {ml_dir}")
        print(f"Library size: {library_size}")

        # ml_table = find_ml_table(outdir_path)
        r_script += f"""
        ml_probability(path_to_table = "{ml_dir}/", table = "{ml_table}", library_size = {library_size});
        """

    if args.html:
        if not args.html_type:
            sys.exit("Error: --html_type must be specified when using --html")
        r_script += f"""
        make_html_summary(path_to_tables = "{outdir_path.resolve()}/", type = "{args.html_type}", ml_plots = {str(args.ml).upper()});
        """

    subprocess.run(["Rscript", "-e", r_script], check=True)

    

def main():
    parser = argparse.ArgumentParser(description='Run MiSiPi.RNA functions with specified parameters.')
    parser.add_argument('--roi', required=True, help='Path to the bed file with regions of interest')
    parser.add_argument('--bam_file', required=True, help='Path to the BAM file of sRNA alignments')
    parser.add_argument('--genome', required=True, help='Path to the genome fasta file')
    parser.add_argument('--plot_output', type=bool, default=True, help='Whether to output plots')
    parser.add_argument('--path_to_RNAfold', required=True, help='Path to RNAfold executable')
    parser.add_argument('--path_to_RNAplot', required=True, help='Path to RNAplot executable')
    parser.add_argument('--pi_pal', default='BlYel', choices=['RdYlBl', 'BlYel', 'yelOrRed', 'MagYel', 'Greens'], help='Palette for piRNA heatmap')
    parser.add_argument('--si_pal', default='RdYlBl', choices=['RdYlBl', 'BlYel', 'yelOrRed', 'MagYel', 'Greens'], help='Palette for siRNA heatmap')
    parser.add_argument('--annotate_region', type=bool, default=True, help='Whether to annotate regions')
    parser.add_argument('--weight_reads', default='None', choices=['Top', 'locus_norm', 'None'], help='Weighting method for read counts')
    parser.add_argument('--gtf_file', help='Path to the GTF file')
    parser.add_argument('--write_fastas', type=bool, default=False, help='Whether to write FASTAs')
    parser.add_argument('--out_type', default='pdf', choices=['pdf', 'png'], help='Output file type')

    parser.add_argument('--miRNA', action='store_true', help='Run miRNA function')
    parser.add_argument('--piRNA', action='store_true', help='Run piRNA function')
    parser.add_argument('--siRNA', action='store_true', help='Run siRNA function')
    parser.add_argument('--run_all', action='store_true', help='Run all RNA functions')
    parser.add_argument('--outdir_name', required=True, help='Custom name for the output directory under misipi_output/')

    # Post-processing options
    parser.add_argument('--ml', action='store_true', help='Run ML classification of loci')
    parser.add_argument('--html', action='store_true', help='Generate HTML summary report')
    parser.add_argument('--html_type', choices=['miRNA', 'piRNA', 'siRNA'], help='Data type for HTML summary')

    args = parser.parse_args()
    
    # Step 1: run core misipi_rna function
    run_misipi_main_rna(args)

    # Step 2: post-processing (optional)
    if args.ml or args.html:
        run_postprocessing_ml_html(args)

if __name__ == "__main__":
    main()
