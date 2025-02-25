#!/usr/bin/env python3
import os
import sys 
import pprint
import glob
import argparse
import re
import pandas
import numpy as np
"""
By: Hans Vasquez-Gross
"""

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--trinotate", help="Input Trinotate V3 file from my previous script")
parser.add_argument("-d", "--da", help="DA master out file with padj values. Will select transcript_ids with <= padj value")
parser.add_argument("-p", "--padj", type=float, default=0.05, help="Padj value to select for from DA table. DEFAULT: 0.05")
parser.add_argument("-l", "--log2foldchange", type=float, help="Log2foldchange to select for from DA table. +- values equal to or greater than; DEFAULT: Null (not filtered)")
parser.add_argument("-o", "--output", help="Filtered output file for only DEG based on padj setting")
args = parser.parse_args()

def print_summarytable(df):
    total = df.count().transcript_id
 
    #pfam counts
    pfam = df.count().Pfam
    eggnm_pfam = df.count().EggNM_PFAMs
    pfam_cols = ["Pfam", "EggNM_PFAMs"]

    #go counts
    pfam_go = df.count().gene_ontology_Pfam
    eggnm_go = df.count().EggNM_GOs
    blastp_go = df.count().gene_ontology_BLASTP
    blastx_go = df.count().gene_ontology_BLASTX
    go_cols = ["gene_ontology_Pfam", "EggNM_GOs", "gene_ontology_BLASTP", "gene_ontology_BLASTX"]

    #descriptions
    eggnm_desc = df.count().EggNM_Description
    blastx_desc = df.count().sprot_Top_BLASTX_hit
    blastp_desc = df.count().sprot_Top_BLASTP_hit
    uniprot_desc = df.count().uniprot
    nr_desc = df.count().nr
    desc_cols = ["EggNM_Description", "sprot_Top_BLASTX_hit", "sprot_Top_BLASTP_hit", "uniprot", "nr"]


    #get multi column counts
    descr = df[desc_cols].dropna(thresh=1)
    gos = df[go_cols].dropna(thresh=1)
    pfams = df[pfam_cols].dropna(thresh=1)

    ##Print output summary
    print("Descriptions")
    print("Total\t%d\t%.2f" % (total, total / total * 100))
    print("EggNM Descriptions\t%d\t%.2f" % (eggnm_desc, eggnm_desc / total * 100))
    print("BlastX Descriptions\t%d\t%.2f" % (blastx_desc, blastx_desc / total * 100))
    print("BlastP Descriptions\t%d\t%.2f" % (blastp_desc, blastp_desc / total * 100))
    print("Uniprot Descriptions\t%d\t%.2f" % (uniprot_desc, uniprot_desc / total * 100))
    print("NR Descriptions\t%d\t%.2f\n" % (nr_desc, nr_desc / total * 100))

    print("GOs")
    print("PFam GOs\t%d\t%.2f" % (pfam_go, pfam_go / total * 100))
    print("EggNM GOs\t%d\t%.2f" % (eggnm_go, eggnm_go / total * 100))
    print("BlastP GOs\t%d\t%.2f" % (blastp_go, blastp_go / total * 100))
    print("BlastX GOs\t%d\t%.2f\n" % (blastx_go, blastx_go / total * 100))

    print("PFams")
    print("PFams\t%d\t%.2f" % (pfam, pfam / total * 100))
    print("EggNM PFams\t%d\t%.2f\n\n" % (eggnm_pfam, eggnm_pfam / total * 100))
    
    
    print("Number of transcripts with at least one description\t%d\t%.2f" % (len(descr), len(descr) / total * 100))
    print("Number of transcripts with at least one GO\t%d\t%.2f" % (len(gos), len(gos) / total * 100))
    print("Number of transcripts with at least one PFAM\t%d\t%.2f" % (len(pfams), len(pfams) / total * 100))

##Main Program
if __name__ == "__main__":


    with open(args.trinotate) as in_handle:
        ##FULL DATA
        df = pandas.read_table(args.trinotate, na_values=".")
        df.index = df["transcript_id"]
      
        print("Full data transcript_id summary")
        print_summarytable(df)

        ##Full data by Geneid
        tgroup = df.groupby(df["gene_id"]).agg(lambda x: np.nan if x.isnull().all() else x.dropna())
        print("\n\nFull data gene_id summary")
        print_summarytable(tgroup)
   
        ##DEG section 
        ##padj and log2fc filtering
        da = pandas.read_table(args.da)
        da_padjfilter = da[da['padj'] <= args.padj]['Chr'].tolist() 
        if args.log2foldchange:
            da_fc = da[abs(da['log2FoldChange']) >= args.log2foldchange]['Chr'].tolist() 
            da_padjfilter = list(set(da_fc) & set(da_padjfilter)) ##Find the intersection of the two lists by converting to set 


        if args.log2foldchange:
            print("\n\nDEG transcript_id summary based on PADJ <= %.2f and abs(log2FoldChange) >= %d" % (args.padj, args.log2foldchange))
            print("Total number of transcript_ids = %d" % len(da_padjfilter))
        else:
            print("\n\nDEG transcript_id summary based on PADJ <= %.2f" % args.padj)
            print("Total number of transcript_ids = %d" % len(da_padjfilter))
        
        deg_df = df[df["transcript_id"].isin(da_padjfilter)]
        print_summarytable(deg_df)
        deg_df.to_csv(args.output, sep="\t", index=False, na_rep=".") 

        ## by gene id
        deg_df_tgroup = deg_df.groupby(df["gene_id"]).agg(lambda x: np.nan if x.isnull().all() else x.dropna())
        deg_df_tgroup_list = deg_df_tgroup['gene_id'].tolist()
        if args.log2foldchange:
            print("\n\nDEG gene_id summary based on PADJ <= %.2f and abs(log2FoldChange) >= %d" % (args.padj, args.log2foldchange))
            print("Total number of gene_ids = %d" % len(deg_df_tgroup_list))
        else:
            print("\n\nDEG gene_id summary based on PADJ <= %.2f" % args.padj)
            print("Total number of gene_ids = %d" % len(deg_df_tgroup_list))

        print_summarytable(deg_df_tgroup)

