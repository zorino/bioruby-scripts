#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maxime@deraspe.net
# review:  	
# date:    	12-12-08
# version: 	0.01
# licence:  	

# Extract Journal References

require 'bio'

f = Bio::FlatFile.auto(ARGV[0])

f.each_entry do |gb|
  ref = gb.references
  ref.each do |r|
    puts r.journal
  end
end
