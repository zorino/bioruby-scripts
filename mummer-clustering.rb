#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# review:  	
# date:    	13-05-13
# version: 	0.01
# licence:  	
#
# Dependency : 	nucmer 3.1
# 		bash (ls, lsof, mv, ps, wc ..)
# 		multifasta-manip-class.rb

require_relative '/multifasta-manip-class'
require 'thread'

# TODO : do it with parallel .. way easier thant thread
require 'parallel'


file = ARGV[0]
@ident = ARGV[1]
@covQ = ARGV[2]
@numT = ARGV[3]

@finalSeqs = []
@seqName = "OUTPUT-seq"


usage = "

$ mummer-clustering.rb <file> <% identity> <% coverage> <number of Thread>

"

if ARGV.length < 4
  abort usage
end


# Split multifasta in singlefasta in the created directory
# Return sorted the list of 
def splitFasta file
  @seqName = file.to_s.gsub(".fasta","").gsub(".fa","")

  if Dir.exists? ("./#{@seqName}")
    puts "Directory #{@seqName}/ exists and is not empty"
  else
    puts "spliting fasta file wait a sec..."
    Dir.mkdir("./#{@seqName}",0770)
    Dir.chdir("./#{@seqName}")
    fasta = FastaParser.new("../#{file}")
    fasta.split
    Dir.chdir("../")
  end

  seqString = `ls -Sr #{@seqName}`
  seqList = seqString.to_s.split("\n")

  File.open("#{@seqName}-sorted.lst","w") do |fout|
    fout.write(seqString.to_s)
  end

  return seqList

end


# Look at the MUMmer result to see if it's positive
# HEADERS is :
# [S1]    [E1]    [S2]    [E2]    [LEN 1] [LEN 2] [% IDY] [COV R] [COV Q] [TAGS]
def incorporateSeqBaseOnMUMmerResult? fprefix
  File.open(fprefix,"r") do |fread|
    while l = fread.gets
      l.chomp!
      res = l.split("\t")
      puts res[0]
      if res[8].to_i >= @covQ.to_i
        return 0
      end
    end
  end
  return 1
end



# Write Seq Name to List File. Check if the Result File is Open with lsof
# Need lsof in bash.
def writeSeqToList outfile, output
  # Check if file already open and wait until is closed
  while ! %x[lsof -F n].split("\n").grep(/#{outfile}/).empty?
    sleep 2
  end
  # Write output to outfile
  File.open("#{outfile}", "a") do |final|
    final.write("#{output}\n")
  end
end



# Iterate through a list of contigs to discriminated redundants
def loopMUMmer iterDown, list

  incorporate = 1
  iterUp = list.length-1

  while incorporate == 1 and iterDown < iterUp
    puts "#{list[iterDown]}"
    nucmerExec = `nucmer -p #{list[iterDown]}-#{list[iterUp]} ../#{@seqName}/#{list[iterUp]} ../#{@seqName}/#{list[iterDown]}`
    showcoords = `show-coords -H -T -I #{@ident} -c #{list[iterDown]}-#{list[iterUp]}.delta > #{list[iterDown]}-#{list[iterUp]}.coords`
    incorporate = incorporateSeqBaseOnMUMmerResult?("#{list[iterDown]}-#{list[iterUp]}.coords")
    iterUp -= 1
  end

  if incorporate == 1
    puts "================================"
    puts "incorporating #{list[iterDown]}"
    puts "================================"
    writeSeqToList("Final-Sequences-IN.list.txt","#{list[iterDown]}")
    `rm #{list[iterDown]}*`
  else
    puts "================================"
    puts "NOT incorporating #{list[iterDown]}"
    puts "================================"
    writeSeqToList("Final-Discarded-Sequences.list.txt","#{list[iterDown]}")    
    `mv "#{list[iterDown]}-#{list[iterUp]}.coords" Homologe-"#{list[iterDown]}".out`
    `rm #{list[iterDown]}*`
  end

end


# Run in Loop Mummer to find redundant sequences
def runMUMmer list

  if ! Dir.exist? "./#{@seqName}-MUMmer/"
    Dir.mkdir("./#{@seqName}-MUMmer/",0770)
  end

  Dir.chdir("./#{@seqName}-MUMmer/")

  iterDown = 0


  # Iterate through all the sequences list and run mummer on it
  while iterDown < (list.length - 1)

    numProc = `ps -C nucmer | tail -n +2 | wc -l`
    while numProc.to_i >= (@numT.to_i)
      sleep 5
      puts "=================\n #{numProc} \n ======================="
      numProc = `ps -C nucmer | tail -n +2 | wc -l`
    end

    fork { loopMUMmer(iterDown,list) }

    iterDown += 1

  end

  # Wait after all alignments (nucmer) to finish
  while numProc.to_i >= 1
    puts "=================\n #{numProc} \n ======================="
    numProc = `ps -C nucmer | tail -n +2 | wc -l`
  end

  # Add Last Sequence in the list as part of the ones to keep
  writeSeqToList("Final-Sequences-IN.list.txt","#{list[list.length-1]}")

end

seqs = splitFasta file
runMUMmer seqs

puts "\n ================ \n ENDEND \n ================\n"
