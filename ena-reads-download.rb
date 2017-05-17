#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	2016-02-15
# version: 	0.01

File.open(ARGV[0], "r") do |f|
  while l = f.gets
    lA = l.chomp.strip
    Dir.mkdir(lA.split(",").join("_"))
    Dir.chdir(lA.split(",").join("_"))
    puts "#{lA.split(",").join("_")} download.."
    lA.split(",").each do |acc|
      # 10 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR210/004/SRR2100054/
      # 9 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR210/SRR210594/
      # ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR552/ERR552836

      if acc[0..2] == "ERR"
        `wget -q -r -nH --cut-dirs=4 --no-parent ftp://ftp.sra.ebi.ac.uk/vol1/fastq/#{acc[0..5]}/#{acc}/`
      elsif acc[0..2] == "SRR"
        if acc.length == 10
          `wget -q -r -nH --cut-dirs=5 --no-parent ftp://ftp.sra.ebi.ac.uk/vol1/fastq/#{acc[0..5]}/00#{acc[-1]}/#{acc}/`
        elsif acc.length == 9
          `wget -q -r -nH --cut-dirs=4 --no-parent ftp://ftp.sra.ebi.ac.uk/vol1/fastq/#{acc[0..5]}/#{acc}/`
        end
      end

    end
    Dir.chdir("../")
  end
end
