#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	2015-06-17
# version: 	0.01


usage = "blastn-extract-hits.rb <blast8 file> <queries length> <subject lenght>"


if ARGV.length < 1
  abort usage
end

results = {}

query_length = {}
subject_length = {}

out_by_contigs = "#{ARGV[0]}.bycontig.sum.tsv"
out_by_hits = "#{ARGV[0]}.byhit.sum.tsv"

File.open(ARGV[1],"r") do |f|
  while l = f.gets
    lA = l.chomp.split("\t")
    query_length[lA[0]] = lA[1].to_i
  end
end

File.open(ARGV[2],"r") do |f|
  while l = f.gets
    lA = l.chomp.split("\t")
    subject_length[lA[0]] = lA[1].to_i
  end
end


#EM_PHG:HG531805 NODE_5_length_97800_cov_28.8995_ID_9    100.00  96      0       0       5892    5797    34655   34750   2.9e-30 135.5
#EM_PHG:DQ466086 NODE_5_length_97800_cov_28.8995_ID_9    96.19   289     11      0       5673    5385    34659   34947   1.8e-102        375.3

File.open(ARGV[0],"r") do |f|
  while l = f.gets

    lA = l.chomp.split("\t")
    next if lA[3].to_i < 500
    next if lA[2].to_f < 70

    key = lA[1]
    hit = lA[0]
    p_id = lA[2].to_f
    length = lA[3].to_i
    qbeg = lA[6].to_i
    qend = lA[7].to_i
    sbeg = lA[8].to_i
    send = lA[9].to_i
    eval = lA[10]
    score = lA[11].to_f

    # key = lA[0]
    # hit = lA[1]
    # p_id = lA[2].to_f
    # length = lA[3].to_i
    # qbeg = lA[8].to_i
    # qend = lA[9].to_i
    # sbeg = lA[6].to_i
    # send = lA[7].to_i
    # eval = lA[10]
    # score = lA[11].to_f

    if qbeg > qend
      tmp = qbeg
      qbeg = qend
      qend = tmp
    end

    if sbeg > send
      tmp = sbeg
      sbeg = send
      send = tmp
    end

    if results.has_key? key

      if results[key].has_key? hit

        results[key][hit] << {sbeg: sbeg,
                              send: send,
                              qbeg: qbeg,
                              qend: qend,
                              length: length,
                              p_id: p_id,
                              score: score}

      else
        results[key][hit] = [{sbeg: sbeg,
                              send: send,
                              qbeg: qbeg,
                              qend: qend,
                              length: length,
                              p_id: p_id,
                              score: score}]
      end


    else

      results[key] = {}
      results[key][hit] = [{sbeg: sbeg,
                            send: send,
                            qbeg: qbeg,
                            qend: qend,
                            length: length,
                            p_id: p_id,
                            score: score}]

    end

  end
end


results_merged = {}

results.each_key do |contig|
  results_merged[contig] = {}

  results[contig].each do |hit, hsp|

    sorted_hits_by_subject = hsp.sort_by { |x| x[:sbeg] }
    sorted_hits_by_query = hsp.sort_by { |x| x[:qbeg] }

    results_merged[contig][hit] = {s_pos: [], q_pos: [], length: 0}
    sorted_hits_by_subject.each do |sh|
      results_merged[contig][hit][:s_pos] << "#{sh[:sbeg]}..#{sh[:send]}"
      results_merged[contig][hit][:q_pos] << "#{sh[:qbeg]}..#{sh[:qend]}"
      results_merged[contig][hit][:length] += sh[:length]
    end

  end

end

merged_hits = {}

out_by_contigs_file = File.open(out_by_contigs, "w")


puts "#Contig\tHit\tNbOfMatch\tCoverage\tContigStart\tContigEnd\tContigPositions\tQueryPositions"
out_by_contigs_file.write("#Contig\tHit\tNbOfMatch\tCoverage\tContigStart\tContigEnd\tContigPositions\tQueryPositions\n")

results_merged.each do |c,hit|

  best_hits = []

  hit.each_key do |h_key|

    l = hit[h_key][:length].to_i
    s_pos = hit[h_key][:s_pos].join(",")
    q_pos = hit[h_key][:q_pos].join(",")
    pos = s_pos.split(/[,|\.\.]/)
    newbeg = pos[0].to_i
    newend = pos[-1].to_i

    toadd = 1
    toreplace = -1

    best_hits.each_with_index do |bHit,index|

      if (newbeg >= bHit[:sbeg] and newbeg <= bHit[:send]) or
        (newend >= bHit[:sbeg] and newend <= bHit[:send]) or
        (newbeg <= bHit[:sbeg] and newend >= bHit[:send])

        if bHit[:match] >= l
          toadd = 0
        else
          toadd = 0
          toreplace = index
        end

      end

    end

    if toadd != 0
      best_hits << {sbeg: newbeg,
                    send: newend,
                    match: l,
                    hit: h_key,
                    s_pos: s_pos,
                    q_pos: q_pos}

    elsif toreplace != -1
      best_hits[toreplace] = {sbeg: newbeg,
                              send: newend,
                              match: l,
                              hit: h_key,
                              s_pos: s_pos,
                              q_pos: q_pos}
    end

  end

  # puts "##### By CONTIG ######################### "
  best_hits.each do |h2|
    cov = ((h2[:match].to_f/subject_length[c])*100).round(2)
    puts "#{c}\t#{h2[:hit]}\t#{h2[:match]}\t#{cov}#{h2[:sbeg]}\t#{h2[:send]}\t#{h2[:s_pos]}\t#{h2[:q_pos]}"
    out_by_contigs_file.write("#{c}\t#{h2[:hit]}\t#{h2[:match]}\t#{cov}\t#{h2[:sbeg]}\t#{h2[:send]}\t#{h2[:s_pos]}\t#{h2[:q_pos]}\n")
  end

  best_hits.each_with_index do |h2,index|
    if ! merged_hits.has_key? h2[:hit]
      merged_hits[h2[:hit]] = {contigs: [],
                               match: 0,
                               pos: []}
    end
    merged_hits[h2[:hit]][:contigs] << c
    merged_hits[h2[:hit]][:match] += h2[:match]
    merged_hits[h2[:hit]][:pos].push(*h2[:q_pos].split(","))
  end

end

out_by_contigs_file.close

out_by_hits_file = File.open(out_by_hits,"w")

puts "##### By HIT ######################### "
puts "#Hit\tMatch\tHitLength\tHitCoverage\tHitPositions\tContigsHit"
out_by_hits_file.write("#Hit\tMatch\tHitLength\tHitCoverage\tHitPositions\tContigsHit\n")

merged_hits.each do |k,v|

  pos = v[:pos].sort_by { |x| x.split("..")[0].to_i }
  l = query_length[k]
  cov = (v[:match].to_f/l).round(2)

  newpos = []

  pos.each do |p|
    curpos = p.split("..")
    if newpos.empty?
      newpos << p
    else
      lastpos = newpos[-1].split("..")
      if curpos[0].to_i <= lastpos[-1].to_i
        newpos[-1] = "#{lastpos[0]}..#{curpos[-1]}"
      else
        newpos << p
      end
    end
  end

  # puts "#{k}\t#{v[:match]}\t#{l}\t#{cov}\t#{pos.join(',')}\t#{v[:contigs].join(',')}"
  puts "#{k}\t#{v[:match]}\t#{l}\t#{cov}\t#{newpos.join(',')}\t#{v[:contigs].join(',')}"
  out_by_hits_file.write("#{k}\t#{v[:match]}\t#{l}\t#{cov}\t#{newpos.join(',')}\t#{v[:contigs].join(',')}\n")

end

out_by_hits_file.close
