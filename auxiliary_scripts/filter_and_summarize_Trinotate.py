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
parser.add_argument("-l", "--list", help="txt file with list of contigs to keep")
parser.add_argument("-t", "--trinotate", help="Trinotate XLS file")
parser.add_argument("-o", "--out", help="Output filtered Trinotate file")
args = parser.parse_args()

##Main Program
if __name__ == "__main__":

   contigList = list()
   geneList = list()
   with open(args.list) as in_handle:
       for line in in_handle:
          contigList.append(line.strip())
          x = re.search(r"TRINITY_.+g\d+", line)
          gene = x.group()
          geneList.append(gene)
   geneList = list(set(geneList))

   print("Number of transcripts to keep:", len(contigList))
   print("Number of genes to keep:", len(geneList))

   out_file = open(args.out, "w")
   #annotatedContigSet = set()
   annotatedContigList = list()
   annotatedGeneList = list()

   annotatedGOContigList = list()
   annotatedGOGeneList = list()

   annotatedKEGGContigList = list()
   annotatedKEGGGeneList = list()

   with open(args.trinotate) as in_handle:
       out_file.write(in_handle.readline())
       for line in in_handle:
          line = line.strip()
          lineArr = line.split("\t")
          gene_id = lineArr[0]
          transcript_id = lineArr[1]
          dataList = lineArr[7:]
          dataList2 = lineArr[2:4]
          dataList = dataList + dataList2 


          if transcript_id in contigList:
              out_file.write(line + "\n")
              res = True in (ele != "." for ele in dataList) 
              x = re.search(r"TRINITY_.+g\d+", transcript_id)
              gene = x.group()
              if "GO" in line:
                  annotatedGOContigList.append(transcript_id)
                  annotatedGOGeneList.append(gene)

              if "KEGG" in line:
                  annotatedKEGGContigList.append(transcript_id)
                  annotatedKEGGGeneList.append(gene)

              if res:
                  annotatedContigList.append(transcript_id)
                  #print(gene)
                  annotatedGeneList.append(gene)

   annotatedContigSet = set(annotatedContigList)
   annotatedGeneSet = set(annotatedGeneList)

   annotatedGOGeneSet = set(annotatedGOGeneList)
   annotatedGOContigSet = set(annotatedGOContigList)

   annotatedKEGGContigSet = set(annotatedKEGGContigList)
   annotatedKEGGGeneSet = set(annotatedKEGGGeneList)
   print("Number of annotated transcripts: %s (%.2f%%)" % (len(annotatedContigSet), len(annotatedContigSet) / len(contigList) * 100))
   print("Number of annotated genes: %s (%.2f%%)" % (len(annotatedGeneSet), len(annotatedGeneSet) / len(geneList) * 100))

   print("Number of annotated GO transcripts: %s (%.2f%%)" % (len(annotatedGOContigSet), len(annotatedGOContigSet) / len(contigList) * 100))
   print("Number of annotated GO genes: %s (%.2f%%)" % (len(annotatedGOGeneSet), len(annotatedGOGeneSet) / len(geneList) * 100))

   print("Number of annotated KEGG transcripts: %s (%.2f%%)" % (len(annotatedKEGGContigSet), len(annotatedKEGGContigSet) / len(contigList) * 100))
   print("Number of annotated KEGG genes: %s (%.2f%%)" % (len(annotatedKEGGGeneSet), len(annotatedKEGGGeneSet) / len(geneList) * 100))

