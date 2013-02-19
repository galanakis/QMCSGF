#! /usr/bin/env ruby

require 'json'
require 'thread'


class Pool
  attr_reader :size
  def initialize(size)
    @size=size
    @jobs=Queue.new
    @pool=Array.new(@size) do |i|
      Thread.new do
        Thread.current[:id]=i
        catch(:exit) do
          loop do
            job, args = @jobs.pop
            job.call(*args)
          end
        end
      end
    end

    def schedule(*args,&block)
      puts args
      @jobs << [block,args]
    end

    def shutdown
      @size.times do
        schedule { throw :exit }
      end
      @pool.map(&:join)
    end
  end
end 


def execute(fin)

  label=File.basename(fin,".*")
  fout=label+'_out.yaml'
  fprocess=fout+'.incomplete'
  if not File.exist?(fout)
    puts "processing #{fout}"
    %x{corv_boson #{fin} > #{fprocess}}
    File.rename(fprocess,fout)
  else
    puts 'File '+fout+' exists. Skipping.'
  end

end


numcores=4

if $0 == __FILE__
  
  if  ARGV.size<3 || ARGV[0]!='-n' || (ARGV[1] =~ /^[0-9]+$/)==nil
    puts "Usage: #{$0} -n <numcores> File1 <File2> <File3> ..."
    exit
  end

  numcores=ARGV[1].to_i

  ARGV.shift(2)

  ARGV.each { |f|
    if not File.exist?(f) 
      puts "File "+f+" does not exist."
      exit
    end
  }

  p = Pool.new(numcores)

  ARGV.each {|f| 

    p.schedule do
      execute(f)
    end

  }

  at_exit { p.shutdown }
end
