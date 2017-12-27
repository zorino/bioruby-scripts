#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	2017-11-01
# version: 	0.01

require 'bio'
require 'parallel'

SELF_PATH=File.expand_path(File.dirname(__FILE__))

def get_species species_file

  species = []
  File.open(species_file) do |f|
    while l = f.gets
      s, nb_of_genomes = l.chomp!.split("\t")
      species.push([s, nb_of_genomes.to_i])
    end
  end
  species_sorted = species.sort_by{ |k,v| v }.reverse

  return species_sorted

end


def read_matrix matrix_f

  puts "# Reading matrix file"
  matrix_full = {}

  File.open(matrix_f, "r") do |f|
    f.gets                      # skip header
    while l = f.gets
      # in THE cluster
      lA = l.chomp!.split("\t")
      _cluster = lA[0]
      matrix_full[_cluster] = []
      lA[1..-1].each do |x|
        next if x == "-"
        x.split(";").each do |_x|
          matrix_full[_cluster].push(_x)
        end
      end
    end
  end

  return matrix_full

end


def read_fasta fasta_f

  puts "# Reading sequence file"
  fasta_full = {}

  flat_f = Bio::FlatFile.auto(fasta_f)
  flat_f.each_entry do |e|
    fasta_full[e.definition] = e
  end

  return fasta_full

end


def extract_consensus_seq gene, matrix, fasta

  # candidates = {}
  best_gene = gene.sort_by { |x,y,z| y }.reverse[0]
  consensus = Bio::Alignment.new()

  # read_matrix matrix_f if ! defined? matrix
  # read_fasta fasta_f if ! defined? fasta

  candidates = matrix[best_gene[0]]
  candidates.each do |c|
    if fasta[c].seq.length == best_gene[2]
      consensus.add_seq(fasta[c].seq, "#{c}")
    end
  end

  consensus_seq = []
  (0..best_gene[2]-1).each do |i|
    pos = []
    consensus.each_seq do |s|
      pos << s[i]
    end
    freq = pos.inject(Hash.new(0)) { |h,v| h[v] += 1; h }
    pos_nt = pos.max_by { |v| freq[v] }
    consensus_seq << pos_nt
  end

  return consensus_seq.join('')

end

# Will extract a consensus for all RNA sequences with the same product name
# used for RNAs but not for CDS since different protein can have the same name (e.g hypothetical protein)
def extract_rna_consensus_cdhit cdhit_f, species_name, proc

  puts "# Extracting RNA consensus sequences"
  fasta_f = cdhit_f.gsub(".cd-hit", "")
  clstr_f = cdhit_f + ".clstr"
  info_f = clstr_f + ".info.tsv"
  matrix_f = clstr_f + ".matrix.tsv"
  path = File.dirname(fasta_f)

  rnas = {}
  File.open(info_f, "r") do |f|
    f.gets
    while l = f.gets
      cluster_id, occurences, length, name, others = l.chomp!.split("\t")
      if ! rnas.has_key? name
        rnas[name] = []
      end
      rnas[name].push([cluster_id, occurences.split("/")[0].to_i, length.to_i])
    end
  end

  if ! Dir.exists? "#{path}/zz_consensus"
    Dir.mkdir "#{path}/zz_consensus"
  end
  if ! Dir.exists? "#{path}/zz_consensus/rna"
    Dir.mkdir "#{path}/zz_consensus/rna"
  end

  matrix_full = read_matrix matrix_f
  fasta_full = read_fasta fasta_f

  Parallel.map(rnas, in_threads: proc) do |k,v|
    puts "  #{k} #{v[0][3]}"
    consensus_seq = extract_consensus_seq v, matrix_full, fasta_full
    if consensus_seq != ""
      bioseq = Bio::Sequence.new(consensus_seq)
      output_file = k.gsub(" ","_") + ".fasta"
      File.open("#{path}/zz_consensus/rna/#{output_file}", "w") do |f|
        f.write(bioseq.output_fasta("#{species_name} #{k}"))
      end
    end
  end

  return ""

end


def extract_cds_consensus_cdhit cdhit_f, species_name, proc

  puts "# Extracting CDS consensus sequences"
  fasta_f = cdhit_f.gsub(".cd-hit", "")
  clstr_f = cdhit_f + ".clstr"
  info_f = clstr_f + ".info.tsv"
  matrix_f = clstr_f + ".matrix.tsv"
  path = File.dirname(fasta_f)

  cds = {}
  File.open(info_f, "r") do |f|
    f.gets
    while l = f.gets
      cluster_id, occurences, length, name, others = l.chomp!.split("\t")
      if ! cds.has_key? name
        cds[cluster_id] = []
      end
      cds[cluster_id].push([cluster_id, occurences.split("/")[0].to_i, length.to_i, name])
    end
  end

  if ! Dir.exists? "#{path}/zz_consensus"
    Dir.mkdir "#{path}/zz_consensus"
  end

  if ! Dir.exists? "#{path}/zz_consensus/cds"

    Dir.mkdir "#{path}/zz_consensus/cds"

    matrix_full = read_matrix matrix_f
    fasta_full = read_fasta fasta_f

    Parallel.map(cds, in_threads: proc) do |k,v|
      puts "  #{k} #{v[0][3]}"
      consensus_seq = extract_consensus_seq v, matrix_full, fasta_full
      if consensus_seq != ""
        bioseq = Bio::Sequence.new(consensus_seq)
        output_file = k.gsub(" ","_") + ".fasta"
        File.open("#{path}/zz_consensus/cds/#{output_file}", "w") do |f|
          f.write(bioseq.output_fasta("#{species_name} #{v[0][3]}"))
        end
      end
    end

  end

  return ""

end


def choose_best_ref_genome genomes, genome_candidates_list, clusters

  # Looking for best ref genome (for contig ordering)
  puts "# Find best ref genome for contig ordering"
  genome_candidates = genomes.select{ |k,v| genome_candidates_list.include?(k) }
  genome_candidates = genome_candidates.to_a.shuffle.to_h
  puts "  Genome Candidates : #{genome_candidates.keys.join(',')}"

  if genome_candidates_list.length == 1 # 1 ref genome
    reference_genome = genome_candidates_list[0]
    puts "  Result (best) = #{reference_genome}"
    return reference_genome
  end

  origin_consensus = genome_get_origin genome_candidates

  # 5 points for a good start
  genome_scores = []
  genome_indices = []

  # order cds and cluster from the origin
  genome_candidates.each_with_index do |(k,g),i|
    _index = g['all_cluster'].index(origin_consensus)
    if _index == 0
      genome_scores[i] = 10
    elsif _index == nil
      genome_scores[i] = -10
    else
      genome_scores[i] = 0
      g['all_cluster'] = g['all_cluster'][_index..-1] + g['all_cluster'][0.._index-1]
      g['all_cds'] = g['all_cds'][_index..-1] + g['all_cds'][0.._index-1]
    end
    genome_indices[i] = 0
  end


  # go to next cds for each genome
  genome_indices.map!{ |i| i=i+1 }

  genome_done = []
  cluster_done = {}
  cluster_done[origin_consensus] = nil
  stop_criterion = false

  # walking for cluster order
  while ! stop_criterion

    cluster_candidates = []

    genome_candidates.each_with_index do |(k,g),i|

      if ! genome_done.include? k and
         ! cluster_done.has_key? g['all_cluster'][genome_indices[i]]

        # puts g['all_cluster'][genome_indices[i]]
        cluster_candidates.push(g['all_cluster'][genome_indices[i]])

        if genome_indices[i]+1 < g['all_cluster'].length
          genome_indices[i] = genome_indices[i] + 1
        else
          genome_done.push(k)
        end

      end
    end

    freq = cluster_candidates.inject(Hash.new(0)) { |h,v| h[v] += 1; h }
    next_clusters = freq.keys.sort_by { |k,v| freq[v] }.reverse

    consensus_clusters = []
    consensus_clusters_indices = []
    consensus_clusters_points = []

    next_clusters.each do |c|
      conserved = true
      c_point = 0.0
      c_indices = []
      genome_candidates.each_with_index do |(k,g),i|
        _index = g['all_cluster'].index(c)
        if _index == nil
          conserved = false
        else
          div = (_index-genome_indices[i]+2)
          if div == 0.0
            div = 100000
          end
          c_point += (1/div)
          # puts "#{_index} - #{genome_indices[i]} + 2 = #{c_point}"
          # puts c_point
          c_indices.push(_index)
        end
      end
      if conserved
        consensus_clusters.push(c)
        consensus_clusters_points.push(c_point)
        consensus_clusters_indices.push(c_indices)
      end
    end

    if ! consensus_clusters.empty?

      next_consensus_index = consensus_clusters_points.each_with_index.max[1]
      consensus_cluster = consensus_clusters[next_consensus_index]
      consensus_indices = consensus_clusters_indices[next_consensus_index]

      genome_candidates.each_with_index do |(k,g),i|
        if consensus_indices[i] > genome_indices[i] - 1
          genome_indices[i] = consensus_indices[i]
          genome_scores[i] += 1
        end
      end

      cluster_done[consensus_cluster] = nil

    else

      genome_indices.map!{ |i| i=i+1 }

    end

    if genome_done.length == genome_indices.length
      stop_criterion = true
    end

  end

  ref_genome_index = genome_scores.each_with_index.max[1]
  reference_genome = genome_candidates.keys[ref_genome_index]

  puts "  Result (best) = #{reference_genome}"
  return reference_genome

end


# return the origin_consensus, the number_of_contigs and genome_indices (start)
def genome_get_origin genomes

  # sort genomes cds
  genomes.each_key do |g|
    genomes[g]['all_cds'].sort!
    genomes[g]['all_cds'].each do |cds|
      genomes[g]['all_cluster'].push(genomes[g][cds])
    end
  end

  # get genome origin
  origin_counts = []
  genomes.each_key do |g|
    origin_counts << genomes[g]['origin_cluster']
    # puts genomes[g]['origin_cluster']
    # puts genomes[g]['origin_cds']
  end
  freq = origin_counts.inject(Hash.new(0)) { |h,v| h[v] += 1; h }
  origin_consensus = origin_counts.max_by { |v| freq[v] }

  return origin_consensus

end


def parse_mauve_output mauve_dir, genome, genomes

  genomes[genome]['contig_order'] = []

  File.open("#{mauve_dir}/#{genome}_contigs.tab","r") do |f|
    in_contig = false
    while l = f.gets

      if l[0..14]  == "Ordered Contigs"
        in_contig = true
        next
      elsif l.chomp!.length == 0
        in_contig = false
      elsif l[0..3] == "type"
        next
      end

      next if in_contig == false

      l.chomp!
      next if l.nil?

      lA = l.split("\t")

      orientation = 1
      orientation = 0 if lA[3] == "complement"
      genomes[genome]['contig_order'].push({contig: lA[1].split(".")[0], orientation: orientation})

    end
  end
  # p genomes[genome]['contig_order']

end


def get_cluster_order genomes, cluster_order, target_genome

  # set and cut (if neccesary) first contig
  target_contig_iterator = 0
  current_contig = genomes[target_genome]['contig_order'][target_contig_iterator]
  target_contig_cluster = genomes[target_genome]
  while ! target_contig_cluster.has_key? current_contig[:contig]
    target_contig_iterator += 1
    current_contig = genomes[target_genome]['contig_order'][target_contig_iterator]
    target_contig_cluster = genomes[target_genome]
  end
  target_contig_cluster[current_contig[:contig]].reverse! if current_contig[:orientation] == 0
  _index = target_contig_cluster[current_contig[:contig]].index(cluster_order[0])
  target_genome_end_clusters = []
  if _index != 0 and _index != nil
    target_genome_end_clusters = target_contig_cluster[current_contig[:contig]][0..(_index-1)]
    target_contig_cluster[current_contig[:contig]] = target_contig_cluster[current_contig[:contig]][_index..-1]
  end

  target_contig_iterator = -1
  current_contig = nil
  target_done = false
  new_cluster_order = []
  _contig_cluster_iterator = 0
  _contig_clusters = []

  cluster_order.each_with_index do |c,i|

    new_cluster_order << c
    next if target_done

    if current_contig == nil
      target_contig_iterator += 1

      if target_contig_iterator > (genomes[target_genome]['contig_order'].length-1)
        target_done = true
        next
      end

      current_contig = genomes[target_genome]['contig_order'][target_contig_iterator]

      while ! target_contig_cluster.has_key? current_contig[:contig] or
           (target_contig_iterator > genomes[target_genome]['contig_order'].length-1)
        target_contig_iterator += 1
        current_contig = genomes[target_genome]['contig_order'][target_contig_iterator]
      end

      _contig_clusters = target_contig_cluster[current_contig[:contig]]
      _contig_clusters.reverse! if current_contig[:orientation] == 0
      _contig_cluster_iterator = 0

    end

    if ! cluster_order.include? (_contig_clusters[_contig_cluster_iterator])
      new_cluster_order << _contig_clusters[_contig_cluster_iterator]
      _contig_cluster_iterator += 1
    end

    if _contig_cluster_iterator > (_contig_clusters.length-1)
      current_contig = nil
    end

  end

  # finish to add target clusters
  if _contig_cluster_iterator < (_contig_clusters.length-1)
    _contig_clusters[_contig_cluster_iterator..-1].each do |_c|
      new_cluster_order << _c if ! new_cluster_order.include? _c
    end
    target_contig_iterator += 1
  end

  while (target_contig_iterator<genomes[target_genome]['contig_order'].length-1)
    current_contig = genomes[target_genome]['contig_order'][target_contig_iterator]
    if target_contig_cluster.has_key? current_contig[:contig]
      _contig_clusters = target_contig_cluster[current_contig[:contig]]
      _contig_clusters.reverse! if current_contig[:orientation] == 0
      _contig_clusters.each do |_c|
        new_cluster_order << _c if ! new_cluster_order.include? _c
      end
    end
    target_contig_iterator += 1
  end

  # insert the end of the first target contig
  cluster_to_insert = []
  target_genome_end_clusters.each do |_c|
    _index = new_cluster_order.index(_c)
    if _index == nil
      cluster_to_insert << _c
    else
      if ! cluster_to_insert.empty?
        new_cluster_order.insert(_index, *cluster_to_insert)
        cluster_to_insert = []
      end
    end
  end

  return new_cluster_order

end


# Extract CDS can take several CDS with same product name
def scaffolding_for_cds_orders species_dir, species_name, proc

  puts "# Extracting CDS with Consensus Cluster Order"
  cdhit_f = species_dir+"/zz_pangenome/CDS.fasta.cd-hit"
  # fasta_f = cdhit_f.gsub(".cd-hit", "")
  clstr_f = cdhit_f + ".clstr"
  info_f = clstr_f + ".info.tsv"
  matrix_f = clstr_f + ".matrix.tsv"
  path = File.expand_path(species_dir)

  # cluster info
  puts "# Reading Cluster Information"
  clusters = {}
  genome_origins = []
  File.open(info_f, "r") do |f|
    f.gets
    while l = f.gets
      cluster_id, occurences, length, name, others = l.chomp!.split("\t")
      if name.include? "DnaA" and name.include? "chromosomal replication"
        genome_origins.push(cluster_id)
      end
      if ! clusters.has_key? name
        clusters[cluster_id] = []
      end
      clusters[cluster_id].push([name, occurences.split("/")[0].to_i, length.to_i])
    end
  end

  # genomes matrix
  puts "# Reading Cluster Matrix"
  genomes = {}
  genomes_list = []
  matrix = {}
  matrix_header = []

  File.open(matrix_f, "r") do |f|

    # get header genomes
    l = f.gets
    lA = l.chomp!.split("\t")
    matrix_header = lA
    lA[1..-1].each do |g|
      genomes_list.push(g)
      genomes[g] = {}
      genomes[g]['all_cds'] = []
      genomes[g]['all_cluster'] = []
      genomes[g]['all_contig'] = []
      genomes[g]['origin_cluster'] = ""
    end

    # build genomes hash
    while l = f.gets
      lA = l.chomp!.split("\t")
      matrix[lA[0]] = lA
      cluster_id = lA[0]
      lA[1..-1].each_with_index do |g,i|
        last_origin = ""
        g.split(/;(?! )/).each do |hit|
          next if hit == "-"
          sample, genome_id, prot_id, prot_locus, gene, product = hit.split("|")
          genomes[genomes_list[i]][prot_locus + "_-_" + genome_id] = cluster_id
          genomes[genomes_list[i]][genome_id] = [] if ! genomes[genomes_list[i]].has_key? genome_id
          genomes[genomes_list[i]]['all_cds'].push(prot_locus + "_-_" + genome_id)
          if ! genomes[genomes_list[i]]['all_contig'].include? genome_id
            genomes[genomes_list[i]]['all_contig'].push(genome_id)
          end
          # genomes[genomes_list[i]]['all_cluster'].push(cluster_id)
          if genome_origins.include? cluster_id
            if last_origin == ""
              genomes[genomes_list[i]]['origin_cluster'] = cluster_id
              genomes[genomes_list[i]]['origin_cds'] = prot_locus + "_-_" + genome_id
              last_origin = prot_locus
            elsif last_origin > prot_locus
              genomes[genomes_list[i]]['origin_cluster'] = cluster_id
              genomes[genomes_list[i]]['origin_cds'] = prot_locus + "_-_" + genome_id
            end
          end
        end
      end
    end
  end

  # order genome cds based on their annotations
  puts "# Ordering CDS / Cluster based on original annotation"
  genomes.each_with_index do |(k,g),i|
    cds = []
    File.open("#{species_dir}/#{k}.gb.fts.tsv") do |f|
      l = f.gets
      while l=f.gets
        lA = l.chomp!.split("\t")
        next if lA[0] != "CDS"
        cds_id = lA[7] + "_-_" + lA[1]
        next if ! g.has_key? cds_id
        cds << cds_id
        g[lA[1]] = [] if ! g.has_key? lA[1]
        g[lA[1]].push(g[cds_id]) # push cluster to contig array
      end
    end
    g['all_cds'].sort_by!{ |x| cds.index(x) }
    g['all_cds'].each do |_c|
      g['all_cluster'] << g[_c]
    end
  end

  ## Chose Best Reference Contig (for contig Ordering)
  # choose genome with the minimal number of contigs
  genomes_nb_contigs = []
  genomes.each_with_index do |(k,g),i|
    genomes_nb_contigs.push(g['all_contig'].length)
  end

  genomes_min_contigs = genomes_nb_contigs.each_index.select{|i| genomes_nb_contigs[i] == genomes_nb_contigs.min}

  ref_genome_candidates = []
  genomes_min_contigs.each do |i|
    ref_genome_candidates.push(genomes_list[i])
  end

  # choose ref genomes for it's conserved order (eliminate genome with rearrangements)
  reference_genome = choose_best_ref_genome genomes, ref_genome_candidates, clusters

  # Mauve contig scaffold based on ref genome
  mauve_dir = "#{species_dir}/zz_pangenome/zz_consensus/zz_mauve_aln"

  if ! Dir.exists? mauve_dir

    puts "# Mauve Aln for Contig Ordering"
    Dir.mkdir(mauve_dir)
    all_mauve_exec = []
    genomes.each_with_index do |(k,g),i|
      if k != reference_genome
        cmd = "java -Xmx500m -cp #{SELF_PATH}/mauve/Mauve.jar org.gel.mauve.contigs.ContigOrderer -output #{mauve_dir}/#{k} -ref #{path}/#{reference_genome}.fasta -draft #{path}/#{k}.fasta"
        all_mauve_exec << {genome: k, cmd: cmd}
      end
    end

    results = Parallel.map(all_mauve_exec, in_threads: proc) do |cmd|

      begin

        sleep rand

        mauve_output = `#{cmd[:cmd]} 2>&1 > /dev/null`
        mauve_genome_dir = "#{mauve_dir}/#{cmd[:genome]}"
        all_aln = Dir.glob("#{mauve_genome_dir}/*")
        aln_array = []
        all_aln.each do |a|
          aln_array << File.basename(a).gsub("alignment","").to_i
        end
        aln_array.sort!.reverse!
        best_aln = "alignment#{aln_array[0]}"
        aln_array[1..-1].each do |a|
          `rm -fr #{mauve_genome_dir}/alignment#{a}`
        end
        `mv #{mauve_genome_dir}/#{best_aln}/*_contigs.tab #{mauve_dir}/`
        `rm -fr #{mauve_genome_dir}`

        puts "  #{cmd[:genome]}"

      rescue
        sleep 2
        retry
      end

    end

  end

  # Mauve alignment parsing
  puts "# Parse Mauve aln output"
  results = Parallel.map(genomes.keys, in_threads: proc) do |g|
    if g != reference_genome
      parse_mauve_output mauve_dir, g, genomes
    end
  end

  puts "# Chromosome walking with contig order"
  # chromosome walking with mauve orientation
  cluster_order = genomes[reference_genome]['all_cluster']
  genomes.each_with_index do |(k,g),i|
    # the contig that start can be split
    next if k == reference_genome
    puts "  #{k}"
    cluster_order = get_cluster_order genomes, cluster_order, k
  end

  puts "# Cluster Order"
  File.open("#{mauve_dir}.tsv", "w") do |fout|
    fout.write(matrix_header.join("\t"))
    fout.write("\n")
    cluster_order.each do |c|
      fout.write(matrix[c].join("\t"))
      fout.write("\n")
    end
  end

  puts "# Number of clusters (#{clusters.length})"
  puts "# Result : zz_mauve_aln.tsv"

  # statistics .. TODO more than that
  File.open("#{mauve_dir}.info.txt", "w") do |fout|
    fout.write("Reference_Genome\t#{reference_genome}")
    fout.write("Number_of_Cluster\t#{clusters.length}")
  end

  return reference_genome

end


usage = "

create-species-model-sum.rb <species name> <species dir> <seq type> <proc> [scaffold]

  seq_type = rRNA, tRNA, CDS
  proc = number of processes

"

if ARGV.length < 3
  abort usage
end

species_name = ARGV[0]
species_dir = ARGV[1]
seq_type = ARGV[2]
proc = 2
if ARGV.length > 3
  proc = ARGV[3].to_i
else
  abort usage
end

scaffold = nil
scaffold = true if ARGV.length > 4

# species = get_species species_file

if seq_type.downcase == "rrna"

  rna_cdhit_file = species_dir+"/zz_pangenome/rRNA.fasta.cd-hit"
  extract_rna_consensus_cdhit rna_cdhit_file, species_name, proc

elsif seq_type.downcase == "trna"

  rna_cdhit_file = species_dir+"/zz_pangenome/tRNA.fasta.cd-hit"
  extract_rna_consensus_cdhit rna_cdhit_file, species_name, proc

elsif seq_type.downcase == "cds"

  cds_hit_file = species_dir+"/zz_pangenome/CDS.fasta.cd-hit"
  extract_cds_consensus_cdhit cds_hit_file, species_name, proc
  if scaffold
    scaffolding_for_cds_orders species_dir, species_name,  proc
  end

end

