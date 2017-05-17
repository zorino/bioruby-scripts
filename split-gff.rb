#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# review:  	
# date:    	13-04-30
# version: 	0.01
# licence:  	

file = ARGV[0] or abort "NO FILE given !! \n NEEDS a multi-gff file like prodigal's output"

gffHead = "## gff-version 3\n"
fout = ""
gffStream = ""

f = File.open(file,"r")

while l = f.gets
  if l[0] == "#"
    if l[1] == "#"
      next
    end

    if fout != ""
      File.open("#{fout}.gff","w") do |fsave|
        fsave.write(gffHead)
        fsave.write(gffStream)
      end
      gffStream.clear
      fout.clear
    end
      gffStream << l
  else
    if fout != ""
      gffStream << l
    else
      fout = l.split("\t")[0]
      gffStream << l
    end
  end
end
