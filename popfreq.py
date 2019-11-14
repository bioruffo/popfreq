#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script builds haplotype frequency tables from population data
obtained from Ensembl. It is useful to determine ancestral and derived alleles.
popfreq is offered under the MIT License.
Please read it here: https://opensource.org/licenses/MIT
Created on Tue May 28 12:46:24 2019
@author: roberto
"""

import pandas as pd

class Globdata:
    '''
    This class aggregates population data obtained from multiple files
    and performs the necessary calculations and visualizations.
    
    '''
    def __init__(self, popdatalist):
        '''
        Load the initial population data.
        :param popdatalist: A list of Popdata objects.
        '''
        assert all(type(x) == Popdata for x in popdatalist)
        self.alldata = popdatalist
        self.rsnames = [item.rsname for item in self.alldata]
        error = self.aggregate()
        assert not error

        
    def aggregate(self):
        '''
        Join all data present in `self.alldata` into a Pandas DataFrame.
        '''
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

    def add(self, popdataitem):
        '''
        Add a Popdata item to the main DataFrame.
        :param popdataitem: a Popdata object.
        '''
        assert type(popdataitem) == Popdata
        assert self.alldata[0].df.index == item.df.index
        self.alldata.append(item)
        self.rsnames.append(item.rsname)
        self.df[item.rsname] = item.rs_series


    def count2(self, rs1, rs2, mask=None):
        '''
        Tabulate genotype counts for each haplotype.
        :param rs1: string containing the rs id of the first SNP.
        :param rs2: string containing the rs id of the second SNP.
        :param mask: mask to be applied on the dataframe before counting.
        '''
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
        
    
    def count3(self, rsmask, rstype, rs1, rs2):
        '''
        Tabulate genotype counts for each haplotype, but only count data where
        the first SNP (rs1) matches a desired genotype.
        :param rsmask: string containing the rs id of the SNP used for masking.
        :param rstype: desired genotype at locus `rsmask`.
        :param rs1: string containing the rs id of the first SNP.
        :param rs2: string containing the rs id of the second SNP.
        
        '''
        print(rsmask+" = "+rstype)
        mask = self.df[rsmask] == rstype
        self.count2(rs1, rs2, mask)

            
    def print_categories(self):
        for rs in self.rsnames:
            print(rs + ': genotypes are ' + ', '.join(set(self.df[rs])))
        print()


class Popdata:
    '''
    Load population data for a specific SNP from an Ensembl-generated file.
    '''
    def __init__(self, rsname, prefix = None):
        if not prefix:
            self.prefix = "./373507-SampleGenotypes-Homo_sapiens_Variation_Sample_"
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
    
    # test data
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
    
