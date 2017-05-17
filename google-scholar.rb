#!/usr/bin/env ruby
require 'net/http'
require 'uri'


def escape_title(s)
    s = s.gsub(/(<\/a>)|(&nbsp;)/i,'')
    s = s.gsub(/&hellip;/i,'...')
    return s
end

def get_citation_list(id, start=0, max_num=10)
    str = '/scholar?start=%d&cites=%s'%[start,id]
    enc_uri = URI.escape(str)
    papers = []
    Net::HTTP.start('scholar.google.com') do |http|
        req = Net::HTTP::Get.new(enc_uri)
        s = http.request(req).body
        localizer = s.index('</b> citing <b>')
        s = s[localizer.to_f+20,s.length()-localizer.to_f+20]
        s = s.gsub(/(<b>)|(<\/b>)/i,'')
        start = 0
        1.upto(max_num){
            pos2 = s.index('</h3>', start)
            pos1 = s.rindex('>', pos2-5)
            title =  s[pos1+1, pos2-pos1-1]
            title = escape_title(title)
            #puts title
            papers << title
            start = pos2 + 10
        }
    end
    return papers
end


puts ' *** Googling paper "%s" *** '%ARGV[0]
enc_uri = URI.escape('/scholar?q="'+ARGV[0]+'"&num=1')

Net::HTTP.start('scholar.google.com') do |http|
    req = Net::HTTP::Get.new(enc_uri)
    s = http.request(req).body
    pos1 = s.index('Cited by ')
    if pos1
        pos2 = s.index('</a>', pos1+9)
        citation_num = Integer(s[pos1+9, pos2-pos1-9])
        pos3 = s.rindex('cites', pos1)
        pos4 = s.index('amp', pos3)
        citation_id = s[pos3+6, pos4-pos3-7]
        puts " -- Google Scholar Citation ID: %s"%citation_id
        puts " -- Cited by the following %d papers:"%citation_num
        pages = citation_num / 10
        papers = []
        0.upto(pages-1){ |num|
            papers += get_citation_list(citation_id, num*10)
        }
        num_in_last_page = citation_num - pages*10
        papers += get_citation_list(citation_id, pages*10, num_in_last_page)
        papers.each_with_index {|p,i|
            puts "[#{i+1}] #{p}"
        }
    else
        puts "no citation found!"
    end
end
