#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# review:  	
# date:    	12-11-25
# version: 	0.01
# licence:  	

require 'bio'

class FastaParser               # Class FastaParser

  def initialize multifastafile	        # Hash = gene: line num
    @multifasta = multifastafile
    abort "File doesn't exist!! \n #{@usage}" unless File.exists?(@multifasta)
    lnum = 0
    @flat_fasta = Bio::FlatFile.auto(@multifasta)
    @fastaH = Hash.new()
    f = File.open(@multifasta, "r")
    f.each_line do |l|
      l.chomp!
      lnum += 1
      if l[0,1] == ">"
        gene = l[1..-1]
        @fastaH[gene] = lnum
        #puts "#{gene}\t#{@fastaH[gene].to_s}"
      end
    end
    f.close
  end


  def printGenes
    @fastaH.each_key do |k|
      puts k
    end
  end


  def getGene(name)

    curFasta = ""
    line = []
    regex = name.downcase.to_s
    @fastaH.each_key do |k|
      if k.downcase.to_s =~ /#{regex}/
        line.push(@fastaH[k].to_i)
      end
    end

    line.sort
    out = ""
    i=0
    c=-1
    inSeq = false

    f = File.open(@multifasta, "r")
    f.each_line do |l|
      i += 1
      if l[0] == ">"
        c += 1
        # if l.downcase.include? name.downcase
        # puts l =~ /regex/
        if l.downcase.to_s =~ /#{regex}/
          inSeq = true
          curFasta = "#{curFasta}#{l}"
        else
          inSeq = false
        end
      elsif i == line[c].to_i
        inSeq = true
        curFasta = "#{curFasta}#{l}"
      elsif inSeq
        curFasta = "#{curFasta}#{l}"
      end
    end

    return curFasta

  end


  def split

    if ! Dir.exist? "single-fasta"
      Dir.mkdir("./single-fasta")
    end

    i=0
    file = File.open("single-fasta/#{i}.fasta","w")

    f = File.open(@multifasta, "r")
    f.each_line do |l|
      if l[0] == ">"
        i += 1
        file.close
        fName = l.split(" ")[0].gsub(">","").gsub("<","").gsub("/","-")
        file = File.open("./single-fasta/#{fName}.fasta","a")
        file.write(l.to_s)
      else
        file.write(l)
      end
    end

    file.close()
    File.delete("./single-fasta/0.fasta")
    f.close()

  end


  def splitNumber numOutput

    numberOfSeq = 0
    f = File.open(@multifasta, "r")
    f.each_line do |l|
      if l[0] == ">"
        numberOfSeq += 1
      end
    end

    numSeqByFile = numberOfSeq/numOutput + 1

    name = "Split-"
    i = 0
    prec = 1

    if numOutput > 10
      prec = 2
    elsif numOutput > 100
      prec = 3
    elsif numOutput > 1000
      prec = 4
    end

    fileC = 1
    # item = sprintf("%#{prec}d",fileC)
    fName = "#{name}#{fileC}.fasta"
    file = File.open(fName,"w")

    f = File.open(@multifasta, "r")
    f.each_line do |l|
      if l[0] == ">"
        i += 1

        if i >= numSeqByFile

          file.close
          fileC += 1
          # item = sprintf("%{prec}d",fileC)
          fName = "#{name}#{fileC}"
          # fName = l.split(" ")[0].gsub(">","")
          file = File.open("#{fName}.fasta","w")
          i = 0

        end

        file.write(l.to_s)

      else
        file.write(l)
      end
    end

    file.close()
    f.close()

  end


  def merge nspacer

    seq = ""
    spacer = "n" * nspacer
    # flat = Bio::FlatFile.auto(@multifasta)
    flat = @flat_fasta

    gff = File.new("#{@multifasta}.merged.gff","w")

    pos = 1

    flat.each do |s|
      seq << "#{s.seq}#{spacer}"
      ftname = s.entry_id
      gff.write("#{@multifasta}\t.\tcontig\t#{pos}\t#{pos+s.seq.length}\t.\t+\t.\tID=#{ftname}\n")
      pos += (s.seq.length + nspacer)
    end

    gff.close
    flat.close

    bioseq = Bio::Sequence.new(seq)
    puts bioseq.output_fasta("#{@multifasta} merged with #{nspacer} n")

  end

  def lengthExtract len

    seq = ""
    bpcount = 0
    i=0
    outname = len.gsub(">","gt-").gsub("<","lt-")
    outfile = @multifasta.to_s.split(".")[0]
    file = File.open("#{outfile}-#{outname}.fasta","w")

    f = File.open(@multifasta, "r")
    f.each_line do |l|
      if l[0] == ">" or f.eof?
        if len.include? ">"
          bpNb = len.gsub(">","")
          if bpcount > bpNb.to_i
            file.write(seq)
          end
        elsif len.include? "<"
          bpNb = len.gsub("<","")
          if bpcount < bpNb.to_i
            file.write(seq)
          end
        else
          abort "Symbol > of < not recognized"
        end
        bpcount = 0
        seq = ""
        seq << l
      else
        bpcount = bpcount + l.length - 1
        seq << l
      end
    end
    f.close
    file.close

  end                           # end lengthExtract

  def sortLength order, file

    # flat = Bio::FlatFile.auto(file)
    flat = @flat_fasta

    lengths = []
    index = 0
    flat.each_entry do |entry|
      lengths << { index: index, length: entry.length.to_i, seq: entry }
      index += 1
    end

    lengths.sort_by! { |e| e[:length] }

    if order == "dsc"
      lengths.reverse!
    end

    lengths.each do |k|
      puts "#{k[:seq]}"
    end

  end

  # Fct: Return part of the sequence specify by loc and strand
  def getSeqLoc seqName, seqLoc, strand, prot=nil, header=nil

    # seq = Bio::FastaFormat.new(getGene(seqName))
    seq = ""
    @flat_fasta.each_entry do |e|
      if e.definition.include? seqName
        seq = e
      end
    end
    @flat_fasta.rewind
    abort if seq == ""

    bioSeq = seq.to_biosequence()
    loc = seqLoc.split("..")
    if strand.to_i == -1
      tmp_seq = Bio::Sequence::NA.new(bioSeq.subseq(loc[0].to_i,loc[1].to_i)).reverse_complement
    elsif strand.to_i == 1
      tmp_seq = Bio::Sequence::NA.new(bioSeq.subseq(loc[0].to_i,loc[1].to_i).to_s.downcase)
    else
      abort "Bad Strand : 1 or -1 needed"
    end

    # convert to amino acid
    if prot==1
      tmp_seq = tmp_seq.translate
    end
    sequence = Bio::Sequence.new(tmp_seq)

    output = ""
    if header != nil
      output = sequence.output_fasta("#{header}",60)
    else
      output = sequence.output_fasta("#{bioSeq.accessions[0]}|#{loc[0]}..#{loc[1]}|#{strand}",60)
    end

    return output

  end

  # Fct: Return unique sequence
  def uniqueSequences file

    seqIDs = {}
    printSeq = 1
    File.open(file,"r") do |f|
      while l = f.gets
        if l[0] == ">"
          key = l.split("\n")[0].split(" ")[0]
          if seqIDs.has_key? key
            printSeq = 0
          else
            seqIDs["#{key}"] = 0
            printSeq = 1
            puts l
          end
        elsif printSeq == 1
          puts l
        end
      end
    end

  end


  def formatSequences file

    File.open(file,"r") do |f|
      while l = f.gets
        if l[0] == ">"
          puts(l)
        else
          seq_len = l.length
          i=0
          while i < seq_len
            puts l[i..i+59]
            i += 60
          end
        end
      end
    end


  end


end                             # end of class FastaParser




# Main #
if __FILE__==$0
  # this will only run if the script was the main, not load'd or require'd
  @usage = 
    "
  multifasta-manip.rb [options] <multifasta>

	options :	geneNames
			getSeq <regex>
			split (will split each sequence)
			splitNb <number of file output> 
			merge <number of N's>
			length <length>
			sortLength <asc||dsc> <out>
			getseq-loc <seqName> <beg..end> <strand>
			getseq-unique
			format (will reformat sequence with 60 width)
"

  abort "No input file ! \n #{@usage}" if ARGV[0].nil?

  fasta = FastaParser.new(ARGV[ARGV.length-1])
  opt = ARGV[0]

  case opt.downcase

  when 'genenames'
    fasta.printGenes

  when 'getseq'
    abort "You didn't enter a Gene Name or ID" if ARGV.length <= 2
    name = ARGV[1]
    gene = fasta.getGene(name)
    puts gene

  when 'split'
    fasta.split

  when 'splitnb'
    fasta.splitNumber ARGV[1].to_i

  when 'merge'
    fasta.merge ARGV[1].to_i

  when 'length'
    fasta.lengthExtract ARGV[1]

  when 'sortlength'
    fasta.sortLength ARGV[1], ARGV[2]

  when 'getseq-loc'
    if ARGV.length < 4
      abort "You need to provide sequence name, beginning, end and strand"
    end
    if ARGV[5]
      puts fasta.getSeqLoc ARGV[1], ARGV[2], ARGV[3], ARGV[4]
    else
      puts fasta.getSeqLoc ARGV[1], ARGV[2], ARGV[3]
    end

  when 'getseq-unique'
    fasta.uniqueSequences ARGV[1]

  when 'format'
    fasta.formatSequences ARGV[1]

  else

    puts @usage

  end


end

