#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# review:  	
# date:    	13-08-16
# version: 	0.01
# licence:  	

require 'uri'
require 'http'


usage = "
  usage :  $ itasser-modeling.rb <Sequence>

"

if ARGV[0]
  seq = File.open(ARGV[0],"r").read
  name = seq.split("\n")[0].gsub(">","")
else
  abort usage
end


uri = URI("http://zhanglab.ccmb.med.umich.edu/cgi-bin/itasser_submit.cgi")
Email = "maxime.deraspe.1@ulaval.ca"
Pass = "IT_e4bm"


req = Net::HTTP::Post.new(uri)

req.set_form_data('SEQUENCE' => seq.to_s, 'REPLY-E-MAIL' => Email, 'password' => Pass, 'TARGET-NAME' => name)

res = Net::HTTP.start(uri.hostname, uri.port) do |http|
  http.request(req)
end

case res
when Net::HTTPSuccess, Net::HTTPRedirection
  puts res.body
else
  res.value
end

