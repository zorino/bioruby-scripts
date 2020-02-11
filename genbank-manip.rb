#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# genbank-manip.rb ---
#
# Filename: genbank-manip.rb
# Description:
# Author: Maxime Déraspe
# Maintainer:
# Created: Tue Nov  6 10:40:15 2012 (-0500)
# Version:
# Last-Updated:
#           By:
#     Update #: 0
#

require 'bio'


class GenbankParser

  def initialize f

    @gbfile = File.new(f, "r").read
    @genbanks = Bio::FlatFile.auto(f)

  end

  ## Methods ##

  # Fct: Get dna sequence
  def getDna (cds, seq)
    loc = cds.locations
    sbeg = loc[0].from.to_i
    send = loc[0].to.to_i
    fasta = Bio::Sequence::NA.new(seq.subseq(sbeg,send))
    # position = "#{sbeg}..#{send}"
    if loc[0].strand == -1
      fasta.reverse_complement!
    end
    dna = Bio::Sequence.auto(fasta)
    return dna
  end


  # Fct: Get Features Sequences
  def getFtsNtSequences prefix="", skip_partial=false
    # @gb.features do |ft|
    prefix << "|" if prefix != "" and prefix[-1] != "|"
    output = ""
    @genbanks.each_entry do |gb|
      gb.each_cds do |ft|
        next if ft.feature == "source"
        accession = gb.accession
        ftH = ft.to_hash
        loc = ft.locations
        skip_this = false
        if skip_partial
          loc.each do |l|
            if l.gt or l.lt
              skip_this = true
              break
            end
          end
        end
        next if skip_this
        gene = []
        product = []
        protId = ""
        gene = ftH["gene"] if !ftH["gene"].nil?
        product = ftH["product"] if !ftH["product"].nil?
        protId = ftH["protein_id"][0] if !ftH["protein_id"].nil?
        locustag = ftH["locus_tag"][0] if !ftH["locus_tag"].nil?
        dna = getDna(ft,gb.to_biosequence)
        seqout = dna.output_fasta("#{prefix}#{accession}|#{protId}|#{locustag}|#{gene[0]}|#{product[0]}",60).chomp!
        output += seqout
        output += "\n"
      end
    end
    return output
  end

  # Fct: Get Features Sequences
  def getFtsProtSequences prefix="", skip_partial=false
    prefix << "|" if prefix != "" and prefix[-1] != "|"
    output = ""
    @genbanks.each_entry do |gb|
      gb.each_cds do |ft|
        accession = gb.accession
        ftH = ft.to_hash
        loc = ft.locations
        skip_this = false
        if skip_partial
          loc.each do |l|
            if l.gt or l.lt
              skip_this = true
              break
            end
          end
        end
        next if skip_this
        gene = []
        product = []
        protId = ""
        if ftH.has_key? "pseudo"
          next
        end
        gene = ftH["gene"] if !ftH["gene"].nil?
        product = ftH["product"] if !ftH["product"].nil?
        protId = ftH["protein_id"][0] if !ftH["protein_id"].nil?
        locustag = ftH["locus_tag"][0] if !ftH["locus_tag"].nil?
        if ftH.has_key? "translation"
          pep = ftH["translation"][0] if !ftH["translation"].nil?
        else
          pep = getDna(ft,gb.to_biosequence).translate
        end
        pepBioSeq = Bio::Sequence.auto(pep)
        seqout = pepBioSeq.output_fasta("#{prefix}#{accession}|#{protId}|#{locustag}|#{gene[0]}|#{product[0]}",60).chomp!
        output += seqout
        output += "\n"
      end
    end
    return output
  end

  # Fct: Get All Mobile Elements Fts
  def  getFtsMobileElement
    output = "Accessiont\tSeqBeg\tSeqEnd\tStrand\tMobile_Element_Type\tNote\n"
    @genbanks.each_entry do |gb|
      gb.features.each do |ft|
        if ft.feature == "mobile_element"
          ftH = ft.to_hash
          loc = ft.locations
          seqBeg = loc[0].from.to_s
          seqEnd = loc[0].to.to_s
          strand = loc[0].strand.to_s
          mobile_element_type = ""
          note = ""
          if ftH.has_key? "mobile_element_type"
            mobile_element_type = ftH['mobile_element_type'][0]
          end
          if ftH.has_key? "note"
            note = ftH['note'][0]
          end
          seqout = "#{gb.accession}\t#{seqBeg}\t#{seqEnd}\t#{strand}\t#{mobile_element_type}\t#{note}"
          output += seqout
          output += "\n"
        end
      end
    end
    return output
  end


  # Fct: Get all features location
  def getFtsLoc location, strand, prefix=""
    prefix << "|" if prefix != "" and prefix[-1] != "|"
    output = ""
    loc = location.split("..")
    protId = ""
    @genbanks.rewind
    gbk_entry = @genbanks.next_entry
    gbk_entry.features do |ft|
      next if ft.feature == "gene"
      ftH = ft.to_hash
      ftloc = ft.locations
      if ftloc[0].from == loc[0].to_i && ftloc[0].to == loc[1].to_i
        gene = []
        product = []
        gene = ftH["gene"] if !ftH["gene"].nil?
        product = ftH["product"] if !ftH["product"].nil?
        protId = ftH["protein_id"][0] if !ftH["protein_id"].nil?
        locustag = ftH["locus_tag"][0] if !ftH["locus_tag"].nil?
        location = "c#{location}" if ftloc[0].strand == -1
        dna = getDna(ft, gbk_entry.to_biosequence)
        seqout = dna.output_fasta("#{prefix}#{@accession}|#{protId}|#{locustag}|#{gene[0]}|#{product[0]}",60).chomp!
        output += seqout
        output += "\n"
        break
      end
    end
    return output
  end

  # Fct: Get a feature at a particular locus_tag
  def getFtLocus prefix=""
    prefix << "|" if prefix != "" and prefix[-1] != "|"
    output = ""
    locustag = ARGV[1]
    protId = ""
    @genbanks.each_entry do |gb|
      gb.features do |ft|
        ftH = ft.to_hash
        ftloc = ft.locations
        if ftH.has_key? "locus_tag" and ftH["locus_tag"][0].eql? locustag
          gene = []
          product = []
          gene = ftH["gene"] if !ftH["gene"].nil?
          product = ftH["product"] if !ftH["product"].nil?
          protId = ftH["protein_id"][0] if !ftH["protein_id"].nil?
          locustag = ftH["locus_tag"][0] if !ftH["locus_tag"].nil?
          if ftloc[0].strand == -1
            location = "c#{ftloc[0].from}..#{ftloc[0].to}"
          else
            location = "#{ftloc[0].from}..#{ftloc[0].to}"
          end
          dna = getDna(ft,gb.to_biosequence)
          seqout = dna.output_fasta("#{prefix}#{gb.accession}|#{protId}|#{locustag}|#{gene[0]}|#{product[0]}",60).chomp!
          output += seqout
          output += "\n"
        end
      end
    end
    return output
  end

  # Fct: get ft dna for a particular protein_id
  def getFtProtID prefix=""
    prefix << "|" if prefix != "" and prefix[-1] != "|"
    output = ""
    protein_id = ARGV[1]
    protId = ""
    @genbanks.each_entry do |gb|
      gb.each_cds do |ft|
        ftH = ft.to_hash
        ftloc = ft.locations
        next if ! ftH.has_key? "protein_id"
        proteinIndex = -1
        i = 0
        ftH["protein_id"].each do |x|
          if x.include? protein_id
            proteinIndex = i
          end
          i+=1
        end
        if proteinIndex != -1
          gene = []
          product = []
          gene = ftH["gene"] if !ftH["gene"].nil?
          product = ftH["product"] if !ftH["product"].nil?
          protId = ftH["protein_id"][proteinIndex] if !ftH["protein_id"].nil?
          locustag = ftH["locus_tag"][proteinIndex] if !ftH["locus_tag"].nil?
          if ftloc[0].strand == -1
            location = "complement(#{ftloc[0].from}..#{ftloc[0].to})"
          else
            location = "#{ftloc[0].from}..#{ftloc[0].to}"
          end
          dna = getDna(ft,gb.to_biosequence)
          seqout = dna.output_fasta("#{gb.accession}|#{protId}|#{locustag}|#{gene[0]}|#{product[0]}",60).chomp!
          output += seqout
          output += "\n"
        end
      end
    end
    return output
  end


  # Fct: Return the full dna sequence
  def getSeq
    @genbanks.rewind
    gbk_entry = @genbanks.next_entry
    bioSeq = gbk_entry.to_biosequence
    sequence = Bio::Sequence.new(bioSeq)
    return sequence.output_fasta("#{bioSeq.accessions[0]}",60)
  end


  # Fct: Return part of the sequence specify by loc
  def getSeqLoc

    loc = location.split("..")
    bioSeq = @genbanks.next_entry.to_biosequence
    if strand.to_i == -1
      sequence = Bio::Sequence.new(bioSeq.subseq(loc[0].to_i,loc[1].to_i).reverse_complement)
    elsif strand.to_i == 1
      sequence = Bio::Sequence.new(bioSeq.subseq(loc[0].to_i,loc[1].to_i))
    else
      abort "Bad Strand : 1 or -1 needed"
    end
    return sequence.output_fasta("#{bioSeq.accessions[0]}|#{loc[0]}..#{loc[1]}|#{strand}",60).chomp!
  end


  # Fct: Get feature of according gene or product matching
  #      exactly the kword (insensitive)
  # Need the sequence attached to genbank file
  def getFt prefix=""
    prefix << "|" if prefix != "" and prefix[-1] != "|"
    output = ""
    kword = ARGV[1]
    gb = @genbanks.each_entry[0]
    seq = gb.to_biosequence
    seqoptions = ""
    for c in 2..ARGV.length-1
      seqoptions += "#{ARGV[c]},"
    end
    # look through all features
    gb.each_cds do |ft|
      ftH = ft.to_hash
      loc = ft.locations
      gene = []
      product = []
      protId = ""
      if (!ftH["gene"].nil? && ftH["gene"][0].downcase.include?(kword.downcase)) or
        (!ftH["product"].nil? && ftH["product"][0].downcase.include?(kword.downcase)) or
        (!ftH["locus_tag"].nil? && ftH["locus_tag"][0].downcase.include?(kword.downcase))
        sbeg = loc[0].from.to_i
        send = loc[0].to.to_i
        fasta = Bio::Sequence::NA.new(seq.subseq(sbeg,send))
        position = "#{sbeg}..#{send}"
        if loc[0].strand == -1
          fasta.reverse_complement!
          position = "c#{position}"
        end
        pep = Bio::Sequence.new(fasta.translate)
        gene = ftH["gene"][0] if !ftH["gene"].nil?
        product = ftH["product"][0] if !ftH["product"].nil?
        protId = ftH["protein_id"][0] if !ftH["protein_id"].nil?
        locustag = ftH["locus_tag"][0] if !ftH["locus_tag"].nil?
        if seqoptions.downcase.include?("pep") or seqoptions.downcase.include?("prot")
          seqout = pep.output_fasta("#{prefix}#{@accession}|#{protId}|#{locustag}|#{gene}|#{product}", 60).chomp!
          output += seqout
          output += "\n"
        else
          dna = Bio::Sequence.auto(fasta)
          seqout = dna.output_fasta("#{prefix}#{@accession}|#{protId}|#{locustag}|#{gene}|#{product}",60).chomp!
          output += seqout
          output += "\n"
        end
      end
    end
    return output
  end


  # Fct: Get all features info
  def getFts
    output = "Feature\tAccessiont\tSeqBeg\tSeqEnd\tStrand\tPartial\tPseudo\tProtID\tLocusTag\tGene\tProt-Product\n"
    @genbanks.each_entry do |gb|
      gb.features do |ft|
        next if ft.feature == "gene"
        ftH = ft.to_hash
        loc = ft.locations
        seqBeg = loc[0].from.to_s
        seqEnd = loc[0].to.to_s
        strand = loc[0].strand.to_s
        partial = 0
        if loc[0].gt or loc[0].lt
          partial = 1
        end
        pseudo = 0
        pseudo = 1 if !ftH["pseudo"].nil?
        gene = []
        product = []
        protId = ""
        gene = ftH["gene"] if !ftH["gene"].nil?
        product = ftH["product"] if !ftH["product"].nil?
        protId = ftH["protein_id"].join(",") if !ftH["protein_id"].nil?
        locustag = ftH["locus_tag"][0] if !ftH["locus_tag"].nil?
        output += "#{ft.feature}\t#{gb.accession}\t#{seqBeg}\t#{seqEnd}\t#{strand}\t#{partial}\t#{pseudo}\t#{protId}\t#{locustag}\t#{gene[0]}\t#{product[0]}\n"
      end
    end
    return output
  end

  # Fct: Get RNA (tRNA or rRNA) features info
  def getFtsRNA
    output =  "Accessiont\tSeqBeg\tSeqEnd\tStrand\tPartial\tProtID\tLocusTag\tGene\tProt-Product\n"
    @genbanks.each_entry do |gb|
      gb.features do |ft|
        next if ! ft.feature.to_s.include? "RNA"
        ftH = ft.to_hash
        loc = ft.locations
        seqBeg = loc[0].from.to_s
        seqEnd = loc[0].to.to_s
        strand = loc[0].strand.to_s
        partial = 0
        if loc[0].gt or loc[0].lt
          partial = 1
        end
        gene = []
        product = []
        protId = ""
        gene = ftH["gene"] if !ftH["gene"].nil?
        product = ftH["product"] if !ftH["product"].nil?
        protId = ftH["protein_id"].join(",") if !ftH["protein_id"].nil?
        locustag = ftH["locus_tag"][0] if !ftH["locus_tag"].nil?
        output += "#{gb.accession}\t#{seqBeg}\t#{seqEnd}\t#{strand}\t#{partial}\t#{protId}\t#{locustag}\t#{gene[0]}\t#{product[0]}\n"
      end
    end
    return output
  end


  # Fct: Get CDS features only
  def getFtsGenes
    output = "Accessiont\tSeqBeg\tSeqEnd\tStrand\tPartial\tProtID\tLocusTag\tGene\tProt-Product\n"
    @genbanks.each_entry do |gb|
      gb.features do |ft|
        next if ft.feature != "gene"
        begin
          ftH = ft.to_hash
          loc = ft.locations
          seqBeg = loc[0].from.to_s
          seqEnd = loc[0].to.to_s
          strand = loc[0].strand.to_s
          gene = []
          product = []
          protId = ""
          gene = ftH["gene"] if !ftH["gene"].nil?
          product = ftH["product"] if !ftH["product"].nil?
          protId = ftH["protein_id"].join(",") if !ftH["protein_id"].nil?
          locustag = ftH["locus_tag"][0] if !ftH["locus_tag"].nil?
          output += "#{gb.accession}\t#{seqBeg}\t#{seqEnd}\t#{strand}\t#{protId}\t#{locustag}\t#{gene[0]}\t#{product[0]}\n"
        rescue
          puts "error"
        end
      end
    end
    return output
  end

  # Fct: Get CDS features only
  def getFtsCDS
    output = "Accessiont\tSeqBeg\tSeqEnd\tStrand\tPartial\tProtID\tLocusTag\tGene\tProt-Product\n"
    @genbanks.each_entry do |gb|
      gb.each_cds do |ft|
        ftH = ft.to_hash
        loc = ft.locations
        seqBeg = loc[0].from.to_s
        seqEnd = loc[0].to.to_s
        strand = loc[0].strand.to_s
        partial = 0
        if loc[0].gt or loc[0].lt
          partial = 1
        end
        gene = []
        product = []
        protId = ""
        gene = ftH["gene"] if !ftH["gene"].nil?
        product = ftH["product"] if !ftH["product"].nil?
        protId = ftH["protein_id"].join(",") if !ftH["protein_id"].nil?
        locustag = ftH["locus_tag"][0] if !ftH["locus_tag"].nil?
        output += "#{gb.accession}\t#{seqBeg}\t#{seqEnd}\t#{strand}\t#{partial}\t#{protId}\t#{locustag}\t#{gene[0]}\t#{product[0]}\n"
      end
    end
    return output
  end


  # Fct: GetTaxon out of GenbankFile
  def getTaxon
    classification = ['Kingdom:', 'Phylum:', 'Class:', 'Order:', 'Family:', 'Genus:', 'Species:', 'Organism:']
    i = 0
    gb = @genbanks.next_entry
    gb.classification.each do |taxon|
      puts classification[i] + "\t" + taxon
      i += 1
    end
    puts classification[i] + "\t" + @gb.organism
  end


  # Fct: Create Locus Tag
  def addLocusTag

    gb = @genbanks.next_entry
    prefix = ARGV[1]
    lab_prefix = ARGV[2]
    if ARGV.length == 3
      lab_prefix = "MyLabName"
    end
    length = (gb.each_cds{ |cds| cds }.length).to_s.length
    ftIndex = 0

    if ARGV.length > 3
      ftIndex = ARGV[3].to_i
    end

    last_locus = nil

    gb.features do |ft|

      next if ft.feature == "source" or ft.feature.include? "misc_"

      index_to_rm = []
      ft.qualifiers.each_with_index do |q,i|
        if q.qualifier == "locus_tag" or q.qualifier == "protein_id"
          index_to_rm << i
        end
      end

      index_to_rm.sort.reverse.each do |i|
        ft.qualifiers.delete_at(i)
      end

      if ft.feature == "gene"
        locusNb = format("%0#{length}d", ftIndex+1)
        newQf = Bio::Feature::Qualifier.new('locus_tag', "#{prefix}_#{locusNb}")
        ft.qualifiers.unshift(newQf)
        ftIndex += 1
        last_locus = "#{prefix}_#{locusNb}"
      end

      if ft.feature == "CDS" or
        ft.feature.include? "RNA"
        if last_locus == nil
          # locusNb = format("%0#{length.to_i-1}d", ftIndex)
          locusNb = format("%0#{length}d", ftIndex)
          newQf = Bio::Feature::Qualifier.new('locus_tag', "#{prefix}_#{locusNb}")
          ft.qualifiers.unshift(newQf)
          if ft.feature == "CDS"
            newQf_protein_id = Bio::Feature::Qualifier.new('protein_id', "gnl|#{lab_prefix}|#{prefix}_#{locusNb}")
            ft.qualifiers.push(newQf_protein_id)
          end
          ftIndex += 1
        else
          newQf_locus_tag = Bio::Feature::Qualifier.new('locus_tag', last_locus)
          ft.qualifiers.unshift(newQf_locus_tag)
          if ft.feature == "CDS"
            newQf_protein_id = Bio::Feature::Qualifier.new('protein_id', "gnl|#{lab_prefix}|#{last_locus}")
            ft.qualifiers.push(newQf_protein_id)
          end
        end
      end

    end

    puts gb.to_biosequence.output(:genbank)

  end


  # Fct: Add Gene for all CDS, tRNA, rRNA
  def addGene

    # prefix = ARGV[1]
    # length = @gb.features.length.to_s.length
    gb = @genbanks.next_entry

    ftIndex = -1
    gene_to_add = {}

    genes_done = {}

    gb.features do |ft|

      ftIndex += 1

      if ft.feature == "gene"
        genes_done[ft.position] = 1
      end

      next if ft.feature != "CDS" and ft.feature != "rRNA" and ft.feature != "tRNA" and ft.feature != "ncRNA"

      next if genes_done.has_key? ft.position

      new_feature = Bio::Feature.new("gene", ft.position)

      ft.qualifiers.each do |q|
        if q.qualifier == "locus_tag" or
          q.qualifier == "gene" or
          q.qualifier == "pseudo" or
          q.qualifier == "note"
          new_feature.append(Bio::Feature::Qualifier.new(q.qualifier,q.value))
        end
      end

      gene_to_add[ftIndex] = new_feature

    end

    gene_to_add.sort_by {|k,v| k}.reverse.each do |k,v|
      gb.features.insert(k,v)
    end

    puts gb.to_biosequence.output(:genbank)

  end


  # Fct: Delete Annotation Feature <feature> <qualifier> <regex>
  def deleteFeature

    gb = @genbanks.next_entry
    feature = ARGV[1]
    qualifier = ARGV[2]
    regex = ARGV[3]

    ftIndex = -1

    gb.features do |ft|

      ftIndex += 1

      next if ft.feature != feature

      # new_feature = Bio::Feature.new("gene", ft.position)
      new_qualifiers = []

      ft.qualifiers.each do |q|
        if q.qualifier == qualifier
          if ! q.value.include? regex
            new_qualifiers << Bio::Feature::Qualifier.new(q.qualifier,q.value)
          end
        else
          new_qualifiers << Bio::Feature::Qualifier.new(q.qualifier,q.value)
        end
      end

      ft.qualifiers = new_qualifiers
      # gene_to_add[ftIndex] = new_feature
      # ft = new_feature

    end

    puts gb.to_biosequence.output(:genbank)

  end


  # Fct: Attached fasta sequence to genbank file
  def attachSeq
    fasta = ARGV[1]
    if ! File.exist? fasta
      abort "Fasta file don't exist"
    end

    @gbfile.each_line do |l|
      if l[0..1] != "//"
        puts l
      end
    end

    puts "ORIGIN"

    ntNum = 0

    File.open(fasta,"r") do |rF|
      while l = rF.gets
        if l[0] == ">"
          next
        end
        ntNum += (l.size - 1)
        printf "%9s ", (ntNum - l.size + 2)
        puts l.scan(/.{1,10}/).join(" ")
      end
    end

    puts "//"

  end

  # Fct: Add hypothetical protein for no product CDS
  def addHP

    gb = @genbanks.next_entry
    gb.each_cds do |cds|
      product = 1
      cds.qualifiers.each do |q|
        product = 0 if q.qualifier == "product"
      end
      if product == 1
        newQf = Bio::Feature::Qualifier.new('product', "hypothetical protein")
        cds.append(newQf)
      end
    end

    puts gb.to_biosequence.output(:genbank)

  end



  # Fct: Split MultiGenbank file into single
  def splitMulti

    Dir.mkdir("single-gbk") if ! Dir.exist? "single-gbk"

    f = nil
    buffer = ""

    inside_gbk = false
    @gbfile.each_line do |l|
      if (l[0..4] == "LOCUS" || l[0..9] == "DEFINITION") and inside_gbk == false
        inside_gbk = true
        # contig = l.scan(/contig.*[0-9]*/)[0]
        if l.include? "seqhdr"    # for prodigal
          accession = l.strip.split(";")[2].gsub("\"","").gsub("seqhdr=","")
        else
          accession = l.strip.split(/\s+/)[1]
        end
        f = File.open("single-gbk/#{accession}.gbk","w")
        f.write(l)
      elsif l[0..1] == "//"
        f.write(l)
        f.close
        inside_gbk = false
      else
        f.write(l)
      end
    end

  end


  # Fct: Difference of CDS between 2 genbank files
  def diffCDS

    gb2 = Bio::GenBank.new(File.read(ARGV[1]))

    new_cds = {}
    gb_diff = {}
    gb2_cds = {}

    gb2.each_cds do |cds|
      location = Bio::Locations.new(cds.position)[0]
      annotation = ""
      locus_tag = ""
      cds.qualifiers.each do |q|
        if q.qualifier == "product"
          annotation = q.value
        elsif q.qualifier == "locus_tag"
          locus_tag = q.value
        end
      end
      if location.strand == -1
        new_cds[location.from] = {location: location,
                                  status: "new",
                                  annotation: annotation,
                                  locus_tag: locus_tag}
      else
        new_cds[location.to] = {location: location,
                                status: "new",
                                annotation: annotation,
                                locus_tag: locus_tag}

      end
    end

    gb = @genbanks.next_entry
    gb.each_cds do |cds|
      location = Bio::Locations.new(cds.position)[0]
      annotation = ""
      locus_tag = ""
      cds.qualifiers.each do |q|
        if q.qualifier == "product"
          annotation = q.value
        elsif q.qualifier == "locus_tag"
          locus_tag = q.value
        end
      end
      key = location.to
      key = location.from if location.strand == -1
      old_cds = {}

      if new_cds.has_key? key

        if new_cds[key][:location].from == location.from and
          new_cds[key][:location].to == location.to

          if new_cds[key][:annotation] != annotation # different annotation
            new_cds[key][:status] = "diff_annotation"

            if location.strand == -1
              __start = location.to
              __end = location.from
            else
              __start = location.from
              __end = location.to
            end
            gb_diff[key] = {status: "diff_annotation",
                            s_start: [__start,__start],
                            s_end: __end,
                            s_strand: location.strand,
                            annotation: ["#{annotation}","#{new_cds[key][:annotation]}"],
                            locus_tag: ["#{locus_tag}","#{new_cds[key][:locus_tag]}"]}
          else
            new_cds[key][:status] = "same"
          end

        else                      # different start
          new_cds[key][:status] = "diff_start"
          if location.strand == -1
            __end = location.from
            __start = [location.to, new_cds[key][:location].to]
          else
            __end = location.to
            __start = [location.from, new_cds[key][:location].from]
          end

          diff_annotation = ["#{annotation}","#{annotation}"]
          if new_cds[key][:annotation] != annotation # different annotation
            diff_annotation = ["#{annotation}","#{new_cds[key][:annotation]}"]
          end
          gb_diff[key] = {status: "diff_start",
                          s_start: __start,
                          s_end: __end,
                          s_strand: location.strand,
                          annotation: diff_annotation,
                          locus_tag: ["#{locus_tag}","#{new_cds[key][:locus_tag]}"]}
        end

      else                        # cds not present in new gb
        if location.strand == -1
          __start = location.to
          __end = location.from
        else
          __start = location.from
          __end = location.to
        end
        gb_diff[key] = {status: "diff_deletion",
                        s_start: [__start,""],
                        s_end: __end,
                        s_strand: location.strand,
                        annotation: ["#{annotation}",""],
                        locus_tag: ["#{locus_tag}",""]}
      end

    end

    new_cds.each do |k,v|
      if v[:status] == "new"
        if v[:location].strand == -1
          __start =  v[:location].to
          __end = v[:location].from
        else
          __start =  v[:location].from
          __end = v[:location].to
        end
        gb_diff[k] = {status: "diff_new",
                      s_start:["",__start],
                      s_end: __end,
                      s_strand: v[:location].strand,
                      annotation: ["",v[:annotation]],
                      locus_tag: ["","#{v[:locus_tag]}"]}
      end
    end

    puts "Difference\tLocusTagRef\tStartRef\tLocusTagNew\tStartNew\tEnd\tStrand\tAnnotationRef\tAnnotationNew"
    gb_diff.sort.to_h.each do |k,v|
      __start = v[:locus_tag][0]
      __start += "\t"
      __start += v[:s_start][0].to_s
      __start += "\t"
      __start += v[:locus_tag][1]
      __start += "\t"
      __start += v[:s_start][1].to_s
      __annotation = v[:annotation].join("\t")
      puts "#{v[:status]}\t#{__start}\t#{v[:s_end]}\t#{v[:s_strand]}\t#{__annotation}"
    end


  end


  # Fct: Merge genbank annotation to fasta Sequence
  def mergeGbkSeq


  end

  # Fct: Shift genbank annotation to fasta sequence
  def shiftGbk

  end







end

## Main ##

if __FILE__==$0

  # Help Usage
  usage =
    "Usage : $ genbank-manip.rb <option> <parameters> <gbk file>

    options:
        // GET
        getfts      return all fts info (except: genes)
        getfts-cds	return all CDS features
        getfts-genes	return all Genes features
        getfts-rna	return all RNA (tRNA rRNA) features
        getfts-dna	return all fts dna sequences
        getfts-prot	return all fts translation
        getfts-mobile	return all mobile_element fts
        getft-search	<kword> [pep] nt by default
        getft-loc	<beg..end>
        getft-locus	<locus-tag>
        getft-protid	<protein_id>
        getseq		return seq
        getseq-loc	<beg..end> <strand> (return seq in location)
        gettaxon	return taxonomy
        getpapers	return articles

        // ADD
        addLocusTag	<locus prefix> <submitter prefix> [start_index] will create locus_tag for all features and protein_id for CDS
        addGene		will add gene feature for each CDS, tRNA, rRNA
        addSeq		<fasta> attached fasta to genbank file
        addHP		will add \"hypothetical protein\" annotation on empty CDS

        deleteFeature <feature> <qualifier> <regex>	will delete the qualifier from feature that contains regex
                                                ex. deleteFeature gene note \'Protein homology\'

        // DIFF (compare 2 gbk file)
                diff-cds	<gbk 2>  compare cds from source gbk file

        // Other
        splitMulti	<multi-genbank> will split it into single gbk

"

  ## Init ##
  opt = ARGV[0].to_s.downcase
  abort "File Not Found ! \n\n#{usage}" if ! File.exist?(ARGV[ARGV.length-1].to_s)
  f = ARGV[ARGV.length-1] or abort usage
  # puts @genbanks.dbclass
  # @gb = Bio::GenBank.new(@gbfile)
  # @accession = @gb.accession
  # @org = @gb.organism
  gbk_parser = GenbankParser.new(f)

  # switch case with the options
  case opt

  when "getfts"
    puts gbk_parser.getFts
  when "getfts-cds"
    puts gbk_parser.getFtsCDS
  when "getfts-genes"
    puts gbk_parser.getFtsGenes
  when "getfts-mobile"
    puts gbk_parser.getFtsMobileElement
  when "getfts-rna"
    puts gbk_parser.getFtsRNA
  when "getfts-dna"
    puts gbk_parser.getFtsNtSequences
  when "getfts-prot"
    puts gbk_parser.getFtsProtSequences
  when "getft-loc"
    if ARGV.length == 4
      strand = ARGV[2]
      location = ARGV[1]
      puts gbk_parser.getFtsLoc location, strand
    else
      abort "You need to specify location and strand !"
    end
  when "getft-locus"
    puts gbk_parser.getFtLocus
  when "getseq"
    puts gbk_parser.getSeq
  when "getseq-loc"
    puts gbk_parser.getSeqLoc
  when "getft-protid"
    puts gbk_parser.getFtProtID
  when "getft-search"
    puts gbk_parser.getFt
  when "diff-cds"
    gbk_parser.diffCDS
  when "gettaxon"
    gbk_parser.getTaxon
  when "addseq"
    gbk_parser.attachSeq
  when "addlocustag"
    gbk_parser.addLocusTag
  when "addgene"
    gbk_parser.addGene
  when "deletefeature"
    gbk_parser.deleteFeature
  when "addhp"
    gbk_parser.addHP
  when "splitmulti"
    gbk_parser.splitMulti
  else
    puts usage
  end

end
