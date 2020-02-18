#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	2020-02-13
# version: 	0.01


usage = "

cd-hit-venn.rb <cd-hit.clstr>

"

if ARGV.length < 1
  abort usage
end

combination_counters = {}
samples_in_cluster = []

File.open(ARGV[0]) do |f|
  while l = f.gets
    if l[0] == ">"
      if !samples_in_cluster.empty?
        samples_in_cluster.uniq!
        samples_in_cluster.sort!
        samples_in_cluster.each do |s|
          if combination_counters.has_key? s
            combination_counters[s] += 1
          else
            combination_counters[s] = 1
          end
        end
        if samples_in_cluster.length > 1
          key = samples_in_cluster.join("_")
          if combination_counters.has_key? key
            combination_counters[key] += 1
          else
            combination_counters[key] = 1
          end
        end
      end
      samples_in_cluster = []
    else
      sample = l.split()[2].split("|")[0].gsub(">","")
      samples_in_cluster << sample
    end
  end
end

combination_counters.each do |k,v|
  puts "#{k} #{v}"
end
