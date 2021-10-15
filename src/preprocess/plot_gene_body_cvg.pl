#!/usr/bin/perl
# vim: syntax=perl tabstop=4 expandtab

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my $options = parse_options();
my @files = @ARGV;
my @tokens = ();
open( OFH, ">" . $$options{rfile} ) or die "Error in opening file, $$options{rfile}, $!\n";
foreach my $file( @files ) {
	open( FH, "<$file" ) or die "Error in oepngint file, $file, $!\n";
	my $line = <FH>;
	chomp $line;
	if( $line =~ /^(\S+)/ ) {
		push @tokens, $1;
	}
	print OFH $line, "\n";
	close FH or die "Error in closing the file, $file, $!\n";
}
print OFH 'data_matrix <- matrix(c(' . join( ",", @tokens ) . '), byrow=T, ncol=100)' . "\n";
print OFH 'rowLabel <- c(';
my @labels = ();
foreach my $token( @tokens ) {
	my $cur_token = $token;
	$cur_token =~ s/\.sorted//;
	$cur_token = "\"${cur_token}\"";
	push @labels, $cur_token;
}
print OFH join(",",@labels) . ')' . "\n";
print OFH "\n";
print OFH 'png("' . $$options{png} . '", width = 8, height = 8, unit="in",res=300)' . "\n";
print OFH 'rc <- cm.colors(ncol(data_matrix))' . "\n";
#print OFH 'heatmap(data_matrix, scale=c("none"),keep.dendro=F, labRow = rowLabel ,Colv = NA,Rowv = NA,labCol=NA,col=cm.colors(256),margins = c(6, 8),ColSideColors = rc,cexRow=1,cexCol=1,xlab="Gene body percentile (5\'->3\')", add.expr=x_axis_expr <- axis(side=1,at=c(1,10,20,30,40,50,60,70,80,90,100),labels=c("1","10","20","30","40","50","60","70","80","90","100")))' . "\n";
print OFH 'junk <- dev.off()' . "\n";
print OFH "\n";
print OFH 'png("' . $$options{curves_png} . '", width = 8, height = 8, unit="in",res=300)' . "\n";
print OFH 'x=1:100' . "\n";
print OFH 'icolor = colorRampPalette(c("#7fc97f","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f"))(' . scalar @labels .')' . "\n";
print OFH 'layout(matrix(c(1,1,1,2,1,1,1,2,1,1,1,2), 4, 4, byrow = TRUE))' . "\n";
print OFH 'plot(x,' . $tokens[0] . ',type=\'l\',xlab="Gene body percentile (5\'->3\')", ylab="Coverage",lwd=0.8,col=icolor[1])' . "\n";
foreach my $index( 1 .. scalar @tokens - 1 ) {
	print OFH 'lines(x,' . $tokens[$index] . ',type=\'l\',col=icolor[' . ($index + 1) . '])' . "\n";
}
print OFH 'par(mar=c(1,0,2,1))' . "\n";
print OFH 'plot.new()' . "\n";
print OFH 'legend(0,1,fill=icolor[1:' . scalar @tokens .'],legend=c(' . join( ",", @labels ) . '))' . "\n";
print OFH 'junk <- dev.off()' . "\n";
close OFH or die "Error in closing the file, $$options{rfile}, $!\n";
exit $?;

sub parse_options {
	my $options = {};
	GetOptions( $options, 'rfile|r=s', 'png|p=s', 'curves_png|c=s', 'help|h' );
	unless( $$options{ 'rfile' } or $$options{ 'png' } or $$options{ 'curves_png' } ) {
		print "Usage: $0 <--rfile|r> <--png|p> <--curves_png|-c>\n";
		exit 1;
	}
	return $options;
}
