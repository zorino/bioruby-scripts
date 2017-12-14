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


def extract_consensus_seq rna, matrix_f, fasta_f

  canditate = []
  best_rna = rna.sort_by { |x,y,z| y }.reverse[0]
  index = best_rna[0].gsub("cluster_","").to_i
  iterator = -1
  File.open(matrix_f, "r") do |f|
    f.gets                      # skip header
    while l = f.gets
      iterator += 1
      break if iterator > index
      next if iterator != index
      # in THE cluster
      lA = l.chomp!.split("\t")
      lA.each do |x|
        canditate.push(*x.split(";"))
      end
    end
  end

  consensus = Bio::Alignment.new()
  flat_f = Bio::FlatFile.auto(fasta_f)
  flat_f.each_entry do |e|
    if canditate.include?(e.definition) and e.seq.length == best_rna[2]
      consensus.add_seq(e.seq, e.definition)
    end
  end

  consensus_seq = []
  (0..best_rna[2]-1).each do |i|
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
# used for RNAs but not for CDS since different protein can have the same name (hypothetical protein)
def extract_rna_consensus_cdhit cdhit_f, species_name

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

  rnas.each do |k, v|
    consensus_seq = extract_consensus_seq v, matrix_f, fasta_f
    if consensus_seq != ""
      bioseq = Bio::Sequence.new(consensus_seq)
      # puts bioseq.output_fasta("#{species_name} #{k}")
      output_file = k.gsub(" ","_") + ".fasta"
      File.open("#{path}/zz_consensus/#{output_file}", "w") do |f|
        f.write(bioseq.output_fasta("#{species_name} #{k}"))
      end
    end
  end

  return ""

end


def contigs_walking_next_cluster genomes, genome_indices, cluster_order, consensus_cluster

  genomes_with_consensus = []
  genomes.each_with_index do |(k,g),i|
    if g['all_cluster'][genome_indices[i]] == consensus_cluster
      genomes_with_consensus.push(k)
      if g['all_cluster'].length == (genome_indices[i] + 1)
        genome_indices[i] = 0
      else
        _cds, _contig = g['all_cds'][genome_indices[i]]
        __cds, __contig = g['all_cds'][genome_indices[i]-1]
        if _contig == _contig
          genome_indices[i] = genome_indices[i] + 1
        end
      end
    end
  end

  genomes.each_with_index do |(k,g), i|

    next if genomes_with_consensus.include? k

    # look for wrong contig direction
    if genome_indices[i] > 1
      if g['all_cluster'][genome_indices[i]-2] == consensus_cluster
        g['all_cluster'] = g['all_cluster'][0..genome_indices[i]-2].reverse + g['all_cluster'][genome_indices[i]-1..-1].reverse
        g['all_cds'] = g['all_cds'][0..genome_indices[i]-2].reverse + g['all_cds'][genome_indices[i]-1..-1].reverse
        g['all_cluster'] = g['all_cluster'][1..-1] + g['all_cluster'][0..0]
        g['all_cds'] = g['all_cds'][1..-1] + g['all_cds'][0..0]
        genome_indices[i] = 0
        next
      end
    else
      if g['all_cluster'][-1] == consensus_cluster
        g['all_cluster'] = g['all_cluster'][0..-2].reverse + g['all_cluster'][-1..-1]
        g['all_cds'] = g['all_cds'][0..-2].reverse + g['all_cds'][-1..-1]
        genome_indices[i] = 0
        next
      end
    end

    _cluster = g['all_cluster'][genome_indices[i]]

    if g['all_cluster'].index(consensus_cluster) != nil # next cluster found

      _cds, _contig = g['all_cds'][genome_indices[i]].split("_-_")
      _index = g['all_cluster'].index(consensus_cluster)
      __cds, __contig = g['all_cds'][_index].split("_-_")

      if _contig == __contig # genomic island or genome recombination (like rRNA)

        present_in_consensus = false
        genomes_with_consensus.each do |_g|
          if genomes[_g]['all_cluster'].index(_cluster) != nil
            present_in_consensus = true
            break
          end
        end

        if ! present_in_consensus

          while _cluster != consensus_cluster
            cluster_order.push(_cluster)
            puts "   ..pusing cluster #{_cluster}"
            genome_indices[i] = genome_indices[i] + 1
            _cluster = g['all_cluster'][genome_indices[i]]
          end

        else
          # do nothing because it's better to catch this cluster
          # in the other genome that follows the consensus
        end

      else                # contig break

        _index = g['all_cluster'].index(consensus_cluster)
        g['all_cluster'] = g['all_cluster'][_index..-1] + g['all_cluster'][0.._index-1]
        g['all_cds'] = g['all_cds'][_index..-1] + g['all_cds'][0.._index-1]
        genome_indices[i] = 1

      end

    else                  # next cluster not found
      # nothing to do cluster absent from genome
    end

  end

end


def contigs_walking  genomes, clusters

  origin_consensus, genome_indices, genome_nb_contigs = genome_get_origin genomes

  # consensus order from clusters
  cluster_order = []
  total_nb_clusters = clusters.keys.length
  current_cluster = origin_consensus

  # order genome's array from the start
  # count the number of contigs
  genomes.each_with_index do |(k,g),i|
    genome_nb_contigs.push(g['all_contigs'].length)
    if g['all_cluster'].include? current_cluster
      index = g['all_cluster'].index current_cluster
      if index != 0
        g['all_cluster'] = g['all_cluster'][index..-1] + g['all_cluster'][0..index-1]
        g['all_cds'] = g['all_cds'][index..-1] + g['all_cds'][0..index-1]
      end
    end
    genome_indices[i] = 0
  end


  while cluster_order.length < total_nb_clusters

    cluster_candidates = []

    genomes.each_with_index do |(k,g),i|
      _cds, _contig = g['all_cds'][genome_indices[i]]
      __cds, __contig = g['all_cds'][genome_indices[i]-1]
      if _contig == _contig
        cluster_candidates.push(g['all_cluster'][genome_indices[i]])
      end
    end

    freq = cluster_candidates.inject(Hash.new(0)) { |h,v| h[v] += 1; h }
    freq_array = freq.values
    freq_array.sort!.reverse!

    # p freq_array
    p freq
    p genome_indices

    if freq.length == 1         # unanimous consensus

      puts "Next consensus cluster - unanimous "
      consensus_cluster = cluster_candidates[0]

      genomes.each_with_index do |(k,g),i|
          genome_indices[i] = genome_indices[i] + 1
      end

    else                        # next clusters candidate all different

      puts "Next consensus cluster - different"
      puts " #Cluster Candidate"

      _index = genome_nb_contigs.rindex(genome_nb_contigs.min)
      # _indices = genome_nb_contigs.each_with_index.select {|e, i| e==genome_nb_contigs.min}.map &:last
      consensus_cluster = cluster_candidates[_index]

      puts " ..#{consensus_cluster}"

      contigs_walking_next_cluster genomes, genome_indices, cluster_order, consensus_cluster

    end

    cluster_order.push(consensus_cluster)

  end

  return cluster_order

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

  all_aln = Dir.glob("#{mauve_dir}/*")
  aln_array = []
  all_aln.each do |a|
    aln_array << File.basename(a).gsub("alignment","").to_i
  end
  aln_array.sort!.reverse!
  best_aln = "alignment#{aln_array[0]}"
  aln_array[1..-1].each do |a|
    `rm -fr #{mauve_dir}/alignment#{a}`
  end

  `mv #{mauve_dir}/#{best_aln}/*_contigs.tab #{mauve_dir}/../`
  `rm -fr #{mauve_dir}`
  new_mauve_dir = File.expand_path("..", "#{mauve_dir}")
  genomes[genome]['contig_order'] = []
  File.open("#{new_mauve_dir}/#{genome}_contigs.tab","r") do |f|
    in_contig = false
    while l = f.gets

      if l[0..14]  == "Ordered Contigs"
        in_contig = true
        next
      elsif l.chomp!.length == 0
        in_contig == false
      elsif l[0..3] == "type"
        next
      end

      next if in_contig == false

      lA = l.chomp!.split("\t")
      orientation = 1
      orientation = 0 if lA[3] == "complement"
      genomes[genome]['contig_order'].push({contig: lA[1], orientation: orientation})

    end
  end

  return best_aln

end

# Extract CDS can take several CDS with same product name
def extract_cds_consensus_cdhit species_dir, species_name, proc

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

  File.open(matrix_f, "r") do |f|

    # get header genomes
    l = f.gets
    lA = l.chomp!.split("\t")
    lA[1..-1].each do |g|
      genomes_list.push(g)
      genomes[g] = {}
      genomes[g]['all_cds'] = []
      genomes[g]['all_cluster'] = []
      genomes[g]['all_contigs'] = []
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
          genomes[genomes_list[i]]['all_cds'].push(prot_locus + "_-_" + genome_id)
          if ! genomes[genomes_list[i]]['all_contigs'].include? genome_id
            genomes[genomes_list[i]]['all_contigs'].push(genome_id)
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


  ## Chose Best Reference Contig (for contig Ordering)
  # choose genome with the minimal number of contigs
  genomes_nb_contigs = []
  genomes.each_with_index do |(k,g),i|
    genomes_nb_contigs.push(g['all_contigs'].length)
  end

  genomes_min_contigs = genomes_nb_contigs.each_index.select{|i| genomes_nb_contigs[i] == genomes_nb_contigs.min}

  ref_genome_candidates = []
  genomes_min_contigs.each do |i|
    ref_genome_candidates.push(genomes_list[i])
  end

  # choose ref genomes for it's conserved order (eliminate genome with rearrangements)
  reference_genome = choose_best_ref_genome genomes, ref_genome_candidates, clusters

  # Mauve contig scaffold based on ref genome
  Dir.mkdir("zz-mauve-aln") if ! File.exists? "./zz-mauve-aln"
  all_mauve_exec = []
  genomes.each_with_index do |(k,g),i|
    if k != reference_genome
      cmd = "java -Xmx500m -cp #{SELF_PATH}/mauve/Mauve.jar org.gel.mauve.contigs.ContigOrderer -output zz-mauve-aln/#{k} -ref #{path}/#{reference_genome}.fasta -draft #{path}/#{k}.fasta"
      all_mauve_exec << {genome: k, cmd: cmd}
    end
  end

  puts "# Orienting contig with Mauve"
  results = Parallel.map(all_mauve_exec, in_threads: proc) do |cmd|
    mauve_output = `#{cmd[:cmd]}`
    best_aln = parse_mauve_output "zz-mauve-aln/#{cmd[:genome]}", cmd[:genome], genomes
    puts "  #{cmd[:genome]}"
  end
  print("\n")

  genomes.each_with_index do |(k,g),i|
    # the contig that start can be split
    puts "  # #{k} contig order "
    p g['contig_order']
  end

  return reference_genome

end



usage = "

create-species-model-sum.rb <species name> <species dir> <seq type> <proc>

  seq_type = rRNA, tRNA, CDS

"

if ARGV.length < 3
  abort usage
end

species_name = ARGV[0]
species_dir = ARGV[1]
seq_type = ARGV[2]
proc = 2
proc = ARGV[3].to_i if ARGV.length > 3

# species = get_species species_file

if seq_type.downcase == "rrna"

  rna_cdhit_file = species_dir+"/zz_pangenome/rRNA.fasta.cd-hit"
  extract_rna_consensus_cdhit rna_cdhit_file, species_name

elsif seq_type.downcase == "trna"

  rna_cdhit_file = species_dir+"/zz_pangenome/tRNA.fasta.cd-hit"
  extract_rna_consensus_cdhit rna_cdhit_file, species_name

elsif seq_type.downcase == "cds"

  cluster_order = extract_cds_consensus_cdhit species_dir, species_name,  proc
  # puts "Finish !!"
  # puts cluster_order.length

end

