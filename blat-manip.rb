#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maxime@deraspe.net
# review:  	
# date:    	12-11-29
# version: 	0.02
# licence:  	


require 'bio'

usage = 
"
  blat-manip [opt] <blat output>
	options:	id=    % identity threshold
			covQ=  % coverage query
			covS=  % coverage subject
			help   print help msg

"

flatreport = Bio::FlatFile.auto(ARGV[ARGV.length-1]) or abort usage

id=0
covQ=0
covS=0

for i in 0..ARGV.length-2
  if ARGV[i].include? "id="
    id = ARGV[i].split("=")[1]
    id = id.to_f/100 if id.to_f > 1
  elsif ARGV[i].include? "covQ="
    covQ = ARGV[i].split("=")[1]
    covQ = covQ.to_f/100 if covQ.to_f > 1
  elsif ARGV[i].include? "covS="
    covS = ARGV[i].split("=")[1]
    covS = covS.to_f/100 if covS.to_f > 1
  else
    abort "Options not known ! \n" + usage
  end
end

good = Hash.new()

flatreport.each_entry do |entry|
  entry.each_hit do |h|
    name = h.query_def
    if  ! good.has_key? name \
      and h.percent_identity >= id.to_f \
      and (h.match.to_f/h.query_len) >= covQ.to_f \
      and (h.match.to_f/h.target_len) >= covS.to_f
      covQuery = (h.match.to_f/h.query_len).round(2)
      covSubject = (h.match.to_f/h.target_len).round(2)
      
      good[name] = "#{h.percent_identity}\t#{covQuery}\t#{covSubject}"
    end
  end
end

flatreport.close

good.each do |k,v|
  puts "#{k}\t#{v}"
end

