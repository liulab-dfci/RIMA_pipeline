#!/usr/bin/perl
# vim: syntax=perl tabstop=4 expandtab

#------------------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: Aug, 1, 2015
#------------------------------------

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my $options = parse_options();
my( $info, $basenames ) = get_info( $$options{ 'file' } );
print_info( $info, $basenames );
exit $?;


sub parse_options {
	my $options = {};
	GetOptions( $options, 'file|f=s@', 'help|h' );
	unless( $$options{ 'file' } ) {
		print "Usage: $0 <--file|-f> [--file|-f]\n";
		exit(1);
	}
	return $options;
}

sub get_info {
	my( $file_list ) = @_;
	my $info = {};
	my $basenames = [];
	for my $file( @$file_list ) {
		my $file_base = basename( $file );
		$file_base =~ s/.txt//;
		push @$basenames, $file_base;
		open( FH, "<$file" ) or die "Error in opening the file, $file, $!\n";
		OUTER_WHILE:
		while( my $line = <FH> ) {
			chomp $line;
			if( $line =~ /^Total\s+Assigned\s+Tags\s+(\d+)/ ) {
				my $total_tags = $1;
				#print $total_tags;
				#my $assigned_tags = <FH>;
				my $sep = <FH>;
				my $header = <FH>;
				while( $line = <FH> ) {
					last OUTER_WHILE if $line =~ /^=/;
					chomp $line;
					my( $feature, $total_bases, $tag_count, $tag_kb ) = (undef, undef, undef, undef);
					if( $line =~ /(\S+)\s+(\d+)\s+(\d+)\s+(\d+)/ ) {
						( $feature, $total_bases, $tag_count, $tag_kb ) = ( $1, $2, $3, $4 );
						#print $feature, $total_bases, $tag_count, $tag_kb;
					}
					my $norm_tag_count = sprintf( "%.2f", ( $tag_count / $total_tags ) );
					$$info{ $file_base }{ $feature } = $norm_tag_count;
				}
			}
		}
		close FH or die "Error in closing the file, $file, $!\n";
	}
	return $info, $basenames;
}

sub print_info {
	my( $info, $basenames ) = @_;
	my $header = join( "\t", ( "Feature", @$basenames ) );
	my @features = keys %{ $$info{ $$basenames[0] } };
	print STDOUT $header, "\n";
	foreach my $feature( @features ) {
		my $cur_info = [];
		foreach my $sample( @$basenames ) {
			push @$cur_info, $$info{ $sample }{ $feature };
		} 
		print STDOUT join( "\t", ( $feature, @$cur_info ) ), "\n"; 
	}
}
