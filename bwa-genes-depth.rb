#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maxime@deraspe.net
# review:
# date:    	12-11-23
# version: 	0.01
#
# Dependency:	multifasta-manip.rb (in $PATH)	

# Result files .tsv
# Genes	GeneLenght	CoverageLength	CoverageDepthAverage	PositionsMatch

# Eventually add mismatch from vcf (variants)

require 'bio'
require_relative './multifasta-manip.rb'

usage = 
" 
  bwa-genes-depth.rb <All-genes.depth> <multifasta>

"

@file = ARGV[0] or abort usage
@@fasta = ARGV[1] or abord usage


class SamtoolDepth              # Class ResultParser

  def initialize(depthFile)
    @depthFile=depthFile
  end

  def getGeneNames
    if defined? @resultString
      @resultString.each do |s|
        puts s.split("\t")[0]
      end      
    end
  end

  def getFastaSeq (geneName)
#    out = `multifasta-manip.rb getGene "#{geneName}" "#{@@fasta}"`
    fasta = FastaParser.new()
    out = fasta.getGene(geneName)
  end

  def printResult	# print tabular complete result
    if defined? @resultString
      @resultString.each do |s|
        puts s
      end
    else
      # initiate variable to nil
      curGene = ""
      curPos = 0
      begPos = 0
      lastPos = 0
      posMatch = ""
      numOfMatch = 0
      cov = 0
      @resultString = Array.new
      @resultString.push "Gene\tGeneLengh\tNTaligned\tCovAligned\tCovDepthAvg\tPositions"
      puts "Gene\tGeneLenght\tNTaligned\tCovAligned\tCovDepthAvg\tPositions"
      File.open(@depthFile, "r") do |inf|
        while result = inf.gets
          resArray = result.split("\t")
          curPos = resArray[1].to_i
          begPos = curPos if begPos.eql? 0
          cov += resArray[2].to_i
          if curGene.eql? resArray[0] # continue if it's the same gene
            if curPos != lastPos+1 and lastPos != 0
              # if we see a hole in the alignment save the region
              posMatch += "#{begPos}..#{lastPos};"
              numOfMatch += (lastPos.to_i - begPos.to_i) + 1
              begPos = curPos
            end
            lastPos = curPos
          else                        # initiate another gene
            if curGene != ""          # save the last gene info
              seqString = getFastaSeq("#{curGene}")
              seq = Bio::FastaFormat.new(seqString)
              seqLength = seq.length
              posMatch += "#{begPos}..#{lastPos};"
              numOfMatch += (lastPos.to_i - begPos.to_i) + 1 
              covDepthAvg = sprintf("%.2f", cov.to_f / numOfMatch)
              covAligned = sprintf("%.2f", numOfMatch.to_f / seqLength)
              @resultString.push "#{curGene}\t#{seqLength}\t#{numOfMatch}\t#{covAligned}\t#{covDepthAvg}\t#{posMatch}"
              puts "#{curGene}\t#{seqLength}\t#{numOfMatch}\t#{covAligned}\t#{covDepthAvg}\t#{posMatch}"
              begPos = 0
            end
            curGene = resArray[0]
            posMatch = ""
            numOfMatch = 0
            lastPos = 0
          end
        end
      end
    end
  end



end                             # End of ResultParser Class



## Main ##

res = SamtoolDepth.new(@file)
res.printResult
#res.getGeneNames
