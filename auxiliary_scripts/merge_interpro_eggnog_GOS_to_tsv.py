import pprint
import argparse
import re
#import gzip

"""
By: Hans Vasquez-Gross
"""

parser = argparse.ArgumentParser()
parser.add_argument("-e", "--eggnog", help="eggnog annotations file")
parser.add_argument("-i", "--interpro", help="interpro")
args = parser.parse_args()

##Main Program
if __name__ == "__main__":
    if(args.eggnog):
        eggnog_fh = open(args.eggnog, "r")
    else:
        raise Exception('EGGNOG Annotations file open failed!')

    if(args.interpro):
        interpro_fh = open(args.interpro, "r")
    else:
        raise Exception('interpro file open failed!')


    dataDict = dict()
    dataDictKO = dict()
    with open(args.eggnog) as in_handle:
        header =  in_handle.readline()
        header =  in_handle.readline()
        header =  in_handle.readline()
        header =  in_handle.readline()
        header =  in_handle.readline()
        for line in in_handle:
            if "GO:" in  line:
               goids = line.split("\t")[9].split(",")
               contig = re.sub(r'\.p[0-9]+', "",  line.split("\t")[0])
               #print(contig, goids)
               dataDict[contig] = goids

            #if "KO:" in  line:
            #   koids = line.split("\t")[19].split(",")
            #   contig = re.sub(r'\.p[0-9]+', "",  line.split("\t")[0])
            #   #print(contig, goids)
            #   dataDictKO[contig] = koids


    #pprint.pprint(dataDict)

    with open(args.interpro) as in_handle:
        for line in in_handle:
            if "GO:" in line:
                #print(line)
                goids = line.split("\t")[13].split("|")
                contig = re.sub(r'\.p[0-9]+', "",  line.split("\t")[0])
                if contig in dataDict:
                    eggnog_goids = dataDict[contig]
                    dataDict[contig] = list(set(eggnog_goids + goids))
                else:
                    dataDict[contig] = goids 

    print("contig_name\tGO")
    for contig in dataDict:
        print("%s\t%s" % (contig, ",".join(dataDict[contig])))
