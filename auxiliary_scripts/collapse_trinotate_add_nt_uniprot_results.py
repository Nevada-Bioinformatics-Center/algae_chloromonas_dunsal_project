#!/usr/bin/env python3
import os
import sys 
import pprint
import glob
import argparse
import re
"""
By: Hans Vasquez-Gross
"""

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--trinotate", help="Input Trinotate file")
parser.add_argument("-u", "--uniprot", help="Diamond output file from Uniprot 90 in blast outfmt6 with specific columns (see README)")
parser.add_argument("-n", "--nr", help="Diamond output file from NR database in blast outfmt6 with specific columns (see README)")
args = parser.parse_args()

def escape_characters_gff3(uniprot_hit):
    uniprot_hit = uniprot_hit.replace(";", "%3B")
    uniprot_hit = uniprot_hit.replace("=", "%3D")
    uniprot_hit = uniprot_hit.replace("&", "%26")
    uniprot_hit = uniprot_hit.replace(",", "%2C")
    return(uniprot_hit)

def print_datatable(dataDict, headers):
    header = "\t".join(headers)

    print("transcript_id\t" + header)
    for transcript_id in dataDict:
        linetxt = transcript_id        
        for colname in headers:
            linetxt += "\t%s"  % dataDict[transcript_id][colname]
        print(linetxt)

##Main Program
if __name__ == "__main__":


    #pprint.pprint(dataStrct)
    dataStrctTrin = dict()
    headers = ["gene_id", "EggNM_Description", "sprot_Top_BLASTX_hit", "sprot_Top_BLASTP_hit", "Pfam", "eggnog", "Kegg", "gene_ontology_BLASTX", "gene_ontology_BLASTP", "gene_ontology_Pfam", "EggNM_GOs", "EggNM_seed_ortholog", "EggNM_eggNOG_OGs", "EggNM_COG_category", "EggNM_KEGG_ko", "EggNM_KEGG_Pathway", "EggNM_KEGG_Module", "EggNM_KEGG_Reaction", "EggNM_KEGG_rclass", "EggNM_PFAMs"]
    transcriptorder = list()

    if args.uniprot:
        headers.append("uniprot")

    if args.nr:
        headers.append("nr")

    with open(args.trinotate) as in_handle:
        for line in in_handle:
            line = line.strip()
            if "#" in line:
                continue
            (gene_id, transcript_id, sprot_Top_BLASTX_hit, infernal, prot_id, prot_coords, sprot_Top_BLASTP_hit, Pfam, SignalP, TmHMM, eggnog, Kegg, gene_ontology_BLASTX, gene_ontology_BLASTP, gene_ontology_Pfam, EggNM_seed_ortholog, EggNM_evalue, EggNM_score, EggNM_eggNOG_OGs, EggNM_max_annot_lvl, EggNM_COG_category, EggNM_Description, EggNM_Preferred_name, EggNM_GOs, EggNM_EC, EggNM_KEGG_ko, EggNM_KEGG_Pathway, EggNM_KEGG_Module, EggNM_KEGG_Reaction, EggNM_KEGG_rclass, EggNM_BRITE, EggNM_KEGG_TC, EggNM_CAZy, EggNM_BiGG_Reaction, EggNM_PFAMs, transcript, peptide) = line.split("\t")
            transcriptorder.append(transcript_id)

            match = re.search(r"Full=([\w\s\/\-{}:]+)", sprot_Top_BLASTX_hit)
            if match:
                sprot_Top_BLASTX_hit = match.group(1)

            match = re.search(r"Full=([\w\s\/\-{}:]+)", sprot_Top_BLASTP_hit)
            if match:
                sprot_Top_BLASTP_hit = match.group(1)

            #print(gene_id, transcript_id, EggNM_Description, sprot_Top_BLASTX_hit, sprot_Top_BLASTP_hit, Pfam, eggnog, Kegg, gene_ontology_BLASTX, gene_ontology_BLASTP, gene_ontology_Pfam, EggNM_GOs, EggNM_seed_ortholog, EggNM_eggNOG_OGs, EggNM_COG_category, EggNM_KEGG_ko, EggNM_KEGG_Pathway, EggNM_KEGG_Module, EggNM_KEGG_Reaction, EggNM_KEGG_rclass, EggNM_PFAMs, sep="\t")

            data = (gene_id, EggNM_Description, sprot_Top_BLASTX_hit, sprot_Top_BLASTP_hit, Pfam, eggnog, Kegg, gene_ontology_BLASTX, gene_ontology_BLASTP, gene_ontology_Pfam, EggNM_GOs, EggNM_seed_ortholog, EggNM_eggNOG_OGs, EggNM_COG_category, EggNM_KEGG_ko, EggNM_KEGG_Pathway, EggNM_KEGG_Module, EggNM_KEGG_Reaction, EggNM_KEGG_rclass, EggNM_PFAMs)
            dataStrctTrin.setdefault(transcript_id, dict())
            for colName in headers:
                dataStrctTrin[transcript_id].setdefault(colName, ".")
         
            counter = 0
            for colDat in data:
                if colDat != ".": 
                    colName = headers[counter]
                    dataStrctTrin[transcript_id][colName] = colDat
                counter += 1
            
    transcriptorder = list(dict.fromkeys(transcriptorder))
    #pprint.pprint(transcriptorder)
    #pprint.pprint(dataStrctTrin)


    ##Parse uniprot data and add to dataStrctTrin
    termBlackList = ["charact", "ypothetical"]
    dataStrct = dict()
    if args.uniprot:
        with open(args.uniprot) as in_handle:
            for line in in_handle:
                (qseqid, qlen, sseqid, slen, evalue, bitscore, stitle, qtitle) = line.strip().split("\t")
                qseqid = qseqid.split("|")[0]
                dataStrct.setdefault(qseqid, list())
                dataStrct[qseqid].append(stitle)

    for seqid in dataStrct:
        if dataStrctTrin[seqid]:
            for descript in dataStrct[seqid]:
                if not any(x in descript for x in termBlackList):
                    dataStrctTrin[seqid]["uniprot"] = descript
                else:
                    dataStrctTrin[seqid]["uniprot"] = "."


    ##Parse nr data and add to dataStrctTrin
    dataStrctNr = dict()
    if args.nr:
        with open(args.nr) as in_handle:
            for line in in_handle:
                (qseqid, qlen, sseqid, slen, evalue, bitscore, stitle, qtitle) = line.strip().split("\t")
                qseqid = qseqid.split("|")[0]
                dataStrctNr.setdefault(qseqid, list())
                dataStrctNr[qseqid].append(stitle)

    for seqid in dataStrctNr:
        if dataStrctTrin[seqid]:
            for descript in dataStrctNr[seqid]:
                #if "charact" and "ypothetical" not in descript:
                if not any(x in descript for x in termBlackList):
                    dataStrctTrin[seqid]["nr"] = descript
                else:
                    dataStrctTrin[seqid]["nr"] = "."


    #pprint.pprint(dataStrctTrin)
    print_datatable(dataStrctTrin, headers)
