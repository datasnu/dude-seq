#!/usr/bin/perl -w

use POSIX;
use Getopt::Std;

getopts("p:", \%opts);
my $bar_link_primer = $opts{p};	$bar_link_primer = "" unless ($bar_link_primer);

my $key_seq = "TCAG";
my $flow_order = "TACG";

my $header_seq = "${key_seq}${bar_link_primer}";


my $dat_line = <STDIN>;	# For Quince format
while ($dat_line = <STDIN>)
{
	chomp($dat_line);
	
	if ($dat_line =~ /(\w+)\s+(\d+)\s+(.+)/)
	{
		my $length = $2;
		my @flow_arr = split(/\s+/, $3);
		my $ret_pos = 8;
		
		if ($bar_link_primer ne "")
		{
			$ret_pos = BaseCall($length, @flow_arr);
		}
		
		print ">$1\n$ret_pos $2 $3\n";
	}
}


sub BaseCall
{
	my ($length, @flowgram) = @_;
	my $ret_seq = "";
	
	for (my $i = 0; $i < $length; $i++)
	{
		my $signal = floor($flowgram[$i] + 0.5);
		my $base = substr($flow_order, $i % 4, 1);
		
		for (my $j = 0; $j < $signal; $j++)
		{
			$ret_seq .= $base;
			
			if ($ret_seq =~ /^$header_seq/)
			{
				my $ret_pos = $i + 1;
				return $ret_pos + (4 - ($ret_pos % 4));	# For converting to multiple of four
			}
		}
	}
}
