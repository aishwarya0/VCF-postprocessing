import os
import  pandas as pd
import vcf
import vcfpy
import sys

# Open file, this will read in the header
class FuncotatorOutputReader:
    if len(sys.argv) == 3:
        def __init__(self,reader,f):
            self.out = open("AnnotatedVCF.tsv", 'w')
            self.reader = vcfpy.Reader.from_path(sys.argv[1])
            self.format= ['AD','DP','AF']
            self.f=sys.argv[2]
            self.header = ['#CHROM', 'POS', 'REF','ID','FILTER','ALT','HUGO_SYMBOL','PROTEIN_CHEANGE','Varinat_type','Variant_classification'] + self.reader.header.samples.names
            self.final_header='\t'.join(self.header)
            self.out.write(self.final_header)
            self.out.write("\n")
            for record in self.reader:
                self.l = [record.CHROM, record.POS, record.REF,record.ID,record.FILTER]
                self.l+= [alt.value for alt in record.ALT]
                self.l+= [n.split("|")[0].replace("[","") for n in record.INFO['FUNCOTATION']]
   
                self.l+= [n.split("|")[18] for n in record.INFO['FUNCOTATION']] 
                self.l += [n.split("|")[5] for n in record.INFO['FUNCOTATION']]
                self.l+= [n.split("|")[7] for n in record.INFO['FUNCOTATION']]
            
                if self.f == self.format[0]:
                    self.l+= [call.data.get('AD')  for call in record.calls]
                elif self.f == self.format[1]:
                    self.l+= [call.data.get('DP')  for call in record.calls]
                elif self.f == self.format[2]:
                    self.l+= [call.data.get('AF')  for call in record.calls]
                self.mat='\t'.join(map(str,self.l))
                self.out.write(self.mat)
                self.out.write("\n")
    else:

        print("Usage:python funcotatorOutputParser.py <funcotator_output.vcf> <Filter=(AF/DP/AD)>")
        sys.exit()
            
my_object = FuncotatorOutputReader(sys.argv[1],sys.argv[2])
