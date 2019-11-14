#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 12:46:24 2019

@author: roberto
"""

import pandas as pd

class Globdata:
    
    def __init__(self, popdatalist):
        self.alldata = popdatalist
        self.rsnames = [item.rsname for item in self.alldata]
        error = self.aggregate()
        assert not error

        
    def aggregate(self):
        print("Checking indexes...")
        for itemno in range(1, len(self.alldata)):
            if not all(self.alldata[0].df.index == self.alldata[itemno].df.index):
                print("Indexes do not match!")
                return 1
        print("Aggregating...")
        self.df = self.alldata[0].df[['Population(s)']]
        for item in self.alldata:
            rsname = item.rsname
            self.df[rsname] = item.rs_series
        return 0

    def count2(self, rs1, rs2, mask=None):
        if mask is not None:
            masked = self.df[mask]
        else:
            masked = self.df
        categories = dict()
        for item in [rs1, rs2]:
            categories[item] = sorted(set(masked[item]))
        print(rs1+' '+rs2)
        line = '\t'
        for cat in categories[rs2]:
            line += (cat + '\t')
        print(line)
        for cat1 in categories[rs1]:
            line = cat1
            for cat2 in categories[rs2]:
                mask = (masked[rs1] == cat1) & (masked[rs2] == cat2)
                line += '\t'+str(len(masked[mask]))
            print(line)
        print()
        
    
    def count3(self, rs1, rstype, rs2, rs3):
        print(rs1+" = "+rstype)
        mask = self.df[rs1] == rstype
        self.count2(rs2, rs3, mask)

            
    def print_categories(self):
        for rs in self.rsnames:
            print(rs + ': ' + ', '.join(set(self.df[rs])))
        print()


class Popdata:
    prefix = "./373507-SampleGenotypes-Homo_sapiens_Variation_Sample_"
    
    def __init__(self, rsname):
        self.rsname = rsname
        self.df = pd.read_csv(self.prefix+rsname+".csv", index_col=0, sep=",")
        gencolname = [x for x in self.df.columns if x.lower().startswith("genotype")][0]
        self.df.rename(columns={gencolname:rsname}, inplace=True)
        self.df = self.df[['Population(s)', rsname]]
        self._categorize()
        self.rs_series = self.df[rsname].apply(self._uniform)
        
    
    def __repr__(self):
        return "Popdata({})".format(self.rsname)
    
    def _categorize(self):
        catset = set(self.df[self.rsname])
        self.cats = dict()
        for item in catset:
            if item[::-1] in self.cats:
                self.cats[item] = item[::-1]
            else:
                self.cats[item] = item
    
    def _uniform(self, genotype):
        return self.cats[genotype]
        

if __name__ == '__main__':
    
    rs = ['rs699',
          'rs4762',
          'rs11122576']

    alldata = []
    for item in rs:
        alldata.append(Popdata(item))
    globdata = Globdata(alldata)
    
    #usage
    globdata.print_categories()
    globdata.count2('rs4762', 'rs11122576')
    globdata.count3('rs699', 'G|G', 'rs4762', 'rs11122576')
    
