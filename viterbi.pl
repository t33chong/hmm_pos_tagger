#!/opt/perl/bin/perl

use strict;
use Getopt::Long;
use Pod::Usage;
$| = 1;

# A Markov Model tagger using Viterbi
# Input:  
#		$epFile - emission probabilities (e.g., massive|JJ 0.001)
#		$tpFile - transition probabilities (e.g., IN|VBD  0.25)
#			Note:  be sure to include starting transition probabilities beginning with 
#			       the sequence boundary (e.g., DT|period  0.30), which acts as our start state

my %options;

# defaults, settings
my $epFile    = "e_prob_matrix.txt";	
my $tpFile    = "t_prob_matrix.txt";	
my $seqBound  = "<s>";	
my $unkWord   = "UNK";	
my $verbose   = 0;
my $check     = 0;
my $oldformat = 0;
my $man = 0;
my $help = 0;

GetOptions(
		'help|?' => \$help,
		man => \$man,
		"emission=s" => \$epFile,
		"transition=s" => \$tpFile,
		"boundary=s" => \$seqBound,
		"unk=s" => \$unkWord,
		"verbose" => \$verbose,
		"check" => \$check,
		"oldformat" => \$oldformat)
or pod2usage(1);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

# model parameters
#
# be careful about the indexation of these parameters.
# read the documentation below
#
# emission probabilities, stored as
#    map from word to
#        ( map from tag to log P(word | tag) )
my %A;
# transition probabilities, stored as
#    map from t_i to
#        ( map from t_{i-1} to log P(t_i | t_{i-1}) )
my %B;
# list of all the tags
my @tagset;
# set of all the words, so we can do unknown word handling.
my %words;

# this appears to be tbe best way to "compute" infinity in Perl
my $neginf = -9**9**9;

if ($oldformat) {
	load_params($epFile, $tpFile);
} else {
	load_params_new($epFile, $tpFile);
}
check_params() if $check;

while (<STDIN>) {
	chomp;

	# read in the sentence to process
	my @inTokens;
	foreach my $tok ( split(/\s+/, $_) ) {
		if (not exists $words{$tok}) {
			print "unknown word '$tok'; replacing with $unkWord\n" if $verbose;
			push @inTokens, $unkWord;
		} else { 
			push @inTokens, $tok;
		}
	}

	# Initial probabilities of the Viterbi trellis are set to log(zero).
	# These change as the trellis is populated and navigated.
	# Period is initialized to log(1.0), since this is our start state
	my @viterbi;
	my @bestState;
	my $maxProb;
	my $maxState;
	foreach ( @tagset ) {
		$viterbi[1]{$_} = $neginf;
	}
	$viterbi[1]{$seqBound } = 0.0;

	# An empty element added at the beginning so everything starts at 1!
	unshift(@inTokens, undef);

	# This is the "induction" phase
	# process all input words
	foreach my $i (1 .. $#inTokens){
		print "'$i: $inTokens[$i]'\n" if $verbose;

		# for all possible tags (current token)
		foreach my $j (1 .. $#tagset){
			$maxProb = $neginf;
			$maxState = undef;

			# for all possible tags (previous token)  (preceding + current = bigram)
			foreach my $k (1 .. $#tagset){

			  # current path probality * emission probability * transition probability
			  my $prob = $viterbi[$i]{$tagset[$k]} +
				  getProbEP( $inTokens[$i], $tagset[$j] ) +
				  getProbTP( $tagset[$j], $tagset[$k] );

			  if ($prob > $neginf){
			      print "\t\ttag:  '$tagset[$j]', previous tag: '$tagset[$k]', $prob\n" if $verbose;
			  }

			  # choose maximal probability and associated state, so far
			  if($prob >= $maxProb){
			    $maxProb = $prob;
			    $maxState = $tagset[$k];
			  }
			}

			# store maximal transition probability and associated state
			$viterbi[$i + 1]{$tagset[$j]} = $maxProb;
			$bestState[$i + 1]{$tagset[$j]} = $maxState;
		}
	}


	# Termination and path-readout
	# Choose end state with the maximal probability
	$maxProb = $neginf;
	foreach my $j ( 1 .. $#tagset ){
		my $prob = $viterbi[$#inTokens + 1]{$tagset[$j]};
		if($prob >= $maxProb){
			$maxProb = $prob;
			$maxState = $tagset[$j];
		}
	}
	my @bestPath;
	$bestPath[$#inTokens + 1] = $maxState;

	# Walk back through the trellis and record path
	for (my $j = $#inTokens; $j > 0; $j--){
		$bestPath[$j] = $bestState[$j + 1]{$bestPath[$j + 1]};
	}

	# Don't print the initial period, since it is a dummy, and hence irrelevant
	shift(@inTokens);

	# Print the best path and associated probability
	if ($maxProb > $neginf) {
		print "\nThe best tag sequence for '@inTokens':\n" if $verbose;
		foreach (2 .. $#bestPath) { print "$bestPath[$_] "; }
		print "\n";
		print "Log probability: " if $verbose;
		print "$maxProb\n" if $verbose;
	} else {
		print "\n" if $verbose;
		print "No path found.\n"
	}
}



#subs

# given two files: the emission probability file and the transition
# probability file, populate the A and B parameters, as well as the
# list of tags and set of known words.
sub load_params($$) {
	my ($emissionProbs, $transitionProbs) = @_;

	my %tags;
	open (ep, "<$emissionProbs") || die "Could not open '$emissionProbs'";
	while(<ep>) {
		($_ =~ /([^\|]+)\|([^\s]+)\s+(.+)$/) && validLogprob(${$A{$1}}{$2} = tolog($3)) || warn "The $epFile has the wrong format";
		++$tags{$2};
		++$words{$1};
	}
	close ep;
	open (tp, "<$transitionProbs") || die "Could not open '$transitionProbs'";
	while(<tp>) {
		($_ =~ /([^\|]+)\|([^\s]+)\s+(.+)$/) && validLogprob(${$B{$1}}{$2} = tolog($3)) || warn "The $tpFile has the wrong format";
		++$tags{$1};
		++$tags{$2};
	}
	close tp;
	@tagset = keys %tags;
	unshift(@tagset, undef);
}

# load parameters from the dense matrix format of hw#2

sub load_params_new($$) {
	my ($emissionProbs, $transitionProbs) = @_;

	my %tags;
	open (tp, "<$transitionProbs") || die "Could not open '$transitionProbs'";
	my $header = <tp>;
	$header =~ s/^\s*//;
	$header =~ s/\s*$//;
	my @headerToks = split /\s+/, $header;
	unshift @headerToks, undef;

	while(<tp>) {
		s/^\s*//; s/\s*$//;
		my @toks = split /\s+/, $_;
		for my $i (1 .. $#toks) {
			${$A{$headerToks[$i]}}{$toks[0]} = tolog($toks[$i]);
			++$tags{$headerToks[$i]};
			++$tags{$toks[0]};
		}
	}
	close tp;
	print "loaded $transitionProbs\n" if $verbose;

	open (ep, "<$emissionProbs") || die "Could not open '$emissionProbs'";
	$header = <ep>;
	$header =~ s/^\s*//;
	$header =~ s/\s*$//;
	@headerToks = split /\s+/, $header;
	unshift @headerToks, undef;
	while(<ep>) {
		s/^\s*//; s/\s*$//;
		my @toks = split /\s+/, $_;
		for my $i (1 .. $#toks) {
			${$B{$toks[0]}}{$headerToks[$i]} = tolog($toks[$i]);
			#printf "P(%s|%s) = %f\n", $toks[0], $headerToks[$i], 2 ** getProbEP($toks[0], $headerToks[$i]);
			++$tags{$headerToks[$i]};
			++$words{$toks[0]};
		}
	}
	close ep;
	print "loaded $emissionProbs\n" if $verbose;
	@tagset = keys %tags;
	unshift(@tagset, undef);
}

# make sure things sum to about 1.

sub check_params() {
	foreach my $tag (@tagset) {
		next if not defined $tag;
		my $sum = 0;
		foreach my $nextTag (@tagset) { $sum += 2 ** getProbTP($nextTag, $tag); }
		print STDERR "WARNING: transition probs for $tag do not sum to 1; sum = $sum\n"
			if abs($sum - 1) > 1e-3;
	}
	foreach my $tag (@tagset) {
		next if not defined $tag;
		next if $tag eq $seqBound;
		my $sum = 0;
		foreach my $word (keys %words) { $sum += 2 ** getProbEP($word, $tag); }
		print STDERR "WARNING: emission probs for $tag do not sum to 1; sum = $sum\n"
			if abs($sum - 1) > 1e-3;
	}
}

# convert a number to log

sub tolog($) {
	my ($val) = @_;
	return $neginf if $val eq "-inf";
	return $val if $val <= 0;
	print STDERR "WARNING: some values are positive; converting to log domain\n";
	return log($val) / log(2);
}


# output emission probabilities
sub getProbEP{
	return $neginf if not defined $B{$_[0]};
	return $neginf if not defined ${$B{$_[0]}}{$_[1]};
	return ${$B{$_[0]}}{$_[1]};
}

# output transition probabilities
sub getProbTP{
	return $neginf if not defined $A{$_[0]};
	return $neginf if not defined ${$A{$_[0]}}{$_[1]};
	return ${$A{$_[0]}}{$_[1]};
}

sub isinf { $_[0] == 9**9**9 || $_[0] == -9**9**9 }
sub isnan { ! defined($_[0] <=> 9**9**9) }

sub validLogprob{ !isnan($_[0]) && ($_[0] <= 0) }

__END__

=head1 NAME

viterbi.pl - Simple code for finding the best sequence

=head1 SYNOPSIS

viterbi.pl [options] < input > output

=head1 OPTIONS

=over8

=item B<-verbose>

Print out information about the search space used in Viterbi

=item B<-check>

Make sure that parameters sum to 1 (+/- 1e-3)

=item B<-unk>

Specify the name used for the unknown word; default = UNK

=item B<-boundary>

Specify the name used for the unknown word; default = <s>

=item B<-emission>

Filename of the emission matrix; default = ep.txt

=item B<-transition>

Filename of the transition matrix; default = tp.txt
