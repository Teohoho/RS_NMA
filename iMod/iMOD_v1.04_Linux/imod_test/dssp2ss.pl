#	Extracts Secondary Structure from a DSSP file.
#
#       Mon 30/11/06
#

if ( !($#ARGV == 1) ) 
{
	print "\nUSAGE: $0 [DSSP-file] [output SS-file]\n";
	print "\n\t[DSSP-file] --> DSSP's text file output (version: CMBI version by ElmK / April 1,2000)\n";
	print "\t[output SS-file] --> Output\n";
	print "\nDescription:\n\tOutputs the Secondary Structure of every residue with a simplified text format.\n";
	print "\nSS file format (.ss files): (Two colums simple format: #res  SS)\n(H=helix, E=strand, T=turn, S=bend, B=Bridge, G=3/10 helix, I=Pi helix, U=unknown)\n";
	exit;
}

my $dssp=$ARGV[0]; # input
my $ssfile=$ARGV[1]; # output

my $debug=0;

my $ss=read_text($dssp,"RESIDUE",1); # Reading after keyword only!
printf "Processing %s. (%d atoms readed).\n",$dssp,$#{@$ss}+1;

open(OUT, ">$ssfile") or die "Unable to write $ssfile: $!";

# Loading Secondary Structure (SS)...
my $cont=0;
foreach(@$ss)
{
	my @line=split( "", $_); # to access characters individually
	if($line[13] ne "!") # avoids errors due to segment jumps...
	{
	  if($line[16] eq " ")
	  { printf OUT "%5d %5s\n",$cont,"U"; }
	  else { printf OUT "%5d %5s\n",$cont,$line[16]; }
	  $cont++;
	}
}
close(OUT);
exit;

##############################################################################

# Loads lines after the keyword "$key" on position "$pos"
sub read_text
{ # Only loads the lines after the keyword
	my $file=shift; # input file
	my $key=shift; # keyword to begin reading
	my $pos=shift; # keyword column 
	print "Loading file $file (key: $key  pos: $pos\n" if $debug;
	my @line;
	my @out;
	my $bool=0;
	my $lines=();
        open (FILE,"$file") || die "Cannot open $file\n";
	my $cont;
        while(<FILE>)
	{
		# This loads the whole line
		if($bool==1) { push( @out, $_ ); }
		else
		{
			chomp;
			if( grep ( !/^#/,$_) )
			{
				push @{$lines->[$cont]}, split ( / +/, trim($_) );
#				print "@{$lines->[$cont]}\n" if $debug;
#				print "$lines->[$cont]->[$pos]\n" if $debug;
				if($lines->[$cont]->[$pos] eq $key) { $bool=1; }

				$cont++;
			}
		}
	}
	return(\@out);
}

# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
