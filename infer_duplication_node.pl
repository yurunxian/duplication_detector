#!/usr/bin/env perl
use Bio::TreeIO;
use Modern::Perl;
use IO::All;
use Data::Dumper;
use Array::Utils qw(:all);
use List::MoreUtils qw(:all);
use Getopt::Long qw(GetOptions);

#define the parameters
my $bootstrap_threshold = 50;
my $sp_list;
my $input;
my $output;
GetOptions(
    'b|bootstrap_threshold=s' => \$bootstrap_threshold,
    'i|input=s' => \$input,
    'o|output=s' => \$output,
    'l|list=s' => \$sp_list
);
if ((!$sp_list)or(!$input)or(!$output)) {
	die("Usage: $0 -b <bootstrap_threshold> -i <gene_tree_with_bs> -o <output_dir> -l <species_list_with_phylogeny_order>\n");
}
my $file_name = "";
if ($input =~ /\//) {my @tmp = split /\//,$input; $file_name = $tmp[-1]} else {$file_name = $input}
open (OUT,">>$output/duplication_result.$bootstrap_threshold.txt");
open (GENE,">>$output/duplication_gene.$bootstrap_threshold.txt");
open (IN,"<$sp_list") or die("no sp list");
my @order = ();
my %clade = ();
my %lineage_count = ();
while (<IN>) {
	chomp;
	my @tmp = split /\s+/,$_;
	$lineage_count{$tmp[1]} += 1;
	push @order,$tmp[0];
	$clade{$tmp[0]} = $tmp[1];
}
my $lineage_number = scalar keys %lineage_count;

#read the newick tree
if (! -e "$input") {die($file_name,"\tno file\n")} else{say $file_name}
my $treeio = Bio::TreeIO->new(-format => 'newick',
			      -file => "$input",
			      -internal_node_id => 'bootstrap');
my $tree = $treeio->next_tree;

#reroot the tree at a random leaf tip of the OUTTEST species
my @leaf_all = ();
my %leaf_id = ();
my $root_clade = 0;
foreach my $node ($tree->get_leaf_nodes) {push @leaf_all,$node->id; $leaf_id{$node->id} = $node}
foreach my $clade (@order){
	if (grep {$_ =~ /^$clade/} @leaf_all) {
		my @outgroup = grep {$_ =~ /^$clade/} @leaf_all;
		@outgroup = sort {$a cmp $b} @outgroup;
		$tree->reroot($leaf_id{$outgroup[0]});
#		say "Reroot the tree at ",$outgroup[0];
		$root_clade = $clade{get_sp_name($outgroup[0])};
		last;
	}
}

#get nodes with a >= "$bootstrap_threshold" bootstrap
my @confident_nodes = grep {(($_->bootstrap)and($_->bootstrap > $bootstrap_threshold))} $tree->get_nodes;

#main function
my %duplicated_gene = ();
my @duplication_result = ();
for (my $i = 0; $i < $lineage_number-1; $i++) {$duplication_result[$i] = 0}
foreach my $node (@confident_nodes){
	#find confident nodes with at least four Descendents which include leaves from at least two clades
	my @children = $node->get_all_Descendents;
	my ($clade_number_first,$min_clade_first) = get_clade_number($node);
	if (($min_clade_first <= $root_clade)or($clade_number_first < 2)) {next}
	my @confident_children = intersect(@confident_nodes,@children);
	my @new_confident_children = ();
	for (my $i = 0; $i < scalar @confident_children; $i++) {
		my ($clade_number,$min_clade) = get_clade_number($confident_children[$i]);
		#avoid paraphylic outgroup caused by the reroot processing
		if (($min_clade > $root_clade)and($clade_number >= 2)) {push @new_confident_children,$confident_children[$i]}
	}
	if (scalar @new_confident_children eq 0) {next}

	#find confident child nodes whose leaves share genes from at least one sp with any other leaves
	#find the node whose leaves contain the most ancient lineage
	my @duplication_result_tmp = ();
	for (my $i = 0; $i < $lineage_number-1; $i++) {$duplication_result_tmp[$i] = 0}
	my %duplicate_child_clade = ();
	foreach my $child (@new_confident_children) {
		my @child_Descendents = $child->get_all_Descendents;
		my @other_Descendents = array_minus(@children,@child_Descendents);
		my ($duplicate_clade,$earliest_clade,$shared_sp_number) = duplication(\@child_Descendents,\@other_Descendents);
		if ($shared_sp_number > 0) {
			$duplication_result_tmp[$earliest_clade-1]++;
			push @{$duplicate_child_clade{$earliest_clade-1}},$child;
		}
	}

	#reduce the redundant node
	my $earliest_duplicate_clade = 0;
	for (my $i = 0; $i < $lineage_number-1; $i++) {
		if ($duplication_result_tmp[$i] > 0) {$earliest_duplicate_clade = $i;last;}
	}
	#find internal node
	my @internal_node = ();
	my $judge = 0;
	foreach my $child (@{$duplicate_child_clade{$earliest_duplicate_clade}}) {
		my $path = get_path_between_nodes($node,$child);
		my @confident_internal_node = intersect(@$path,@confident_children);
		if (@confident_internal_node eq 0) {
			$duplication_result[$earliest_duplicate_clade] += 1;
			$judge++;
			last;
		}
	}
	if ($judge > 0) {
		my $duplicated_gene = get_duplicated_gene($node);
		push @{$duplicated_gene{$earliest_duplicate_clade}},@$duplicated_gene;
	}
}

#output
print OUT $file_name;
for (my $i = 1; $i < $lineage_number-1; $i++) {print OUT "\t",$duplication_result[$i]}
print OUT "\n";
print GENE $file_name;
for (my $i = 1; $i < $lineage_number-1; $i++) {print GENE "\t",join(",",(uniq(@{$duplicated_gene{$i}})))}
print GENE "\n";

sub get_sp_name {
	my $gene = shift;
	if ($gene =~ /^(\S+?)\_/) {return $1}
}

sub get_clade_number {
	my $node = shift;
	my @children = $node->get_all_Descendents;
	my %children_leaf = ();
	my $leaf_number = 0;
	for (@children) {
		if ($_->id) {
			$leaf_number++;
			$children_leaf{$clade{get_sp_name($_->id)}}++;
		}
	}
	my @keys = keys %children_leaf;
	my ($clade_number,$min_clade) = (0,0,);
	$clade_number = scalar @keys;
	if ($clade_number > 0) {
		my @tmp = sort {$a <=> $b} keys %children_leaf;
		$min_clade = $tmp[0];
	}
	return ($clade_number,$min_clade);
}

sub duplication {
	my ($tmp1,$tmp2) = @_;
	my @nearest = @$tmp1;
	my @other = @$tmp2;
	my (%nearest_leaf,%other_leaf);
	for (@nearest) {if ($_->id) {$nearest_leaf{get_sp_name($_->id)}++}}
	for (@other) {if ($_->id) {$other_leaf{get_sp_name($_->id)}++}}
	my ($duplicate_clade,$earliest_clade,$shared_sp_number) = (0,0,0);
	for(@order) {if ((exists $nearest_leaf{$_})and(exists $other_leaf{$_})) {$shared_sp_number++}}
	for(@order) {if ((exists $nearest_leaf{$_})and(exists $other_leaf{$_})) {$duplicate_clade = $clade{$_}; last;}}
	for(@order) {if (exists $nearest_leaf{$_}) {$earliest_clade = $clade{$_}; last;}}
	return ($duplicate_clade,$earliest_clade,$shared_sp_number);
}

sub get_duplicated_gene {
	my $node = shift;
	my @tips = $node->get_all_Descendents;
	my @output = ();
	for (@tips) {if ($_->id) {push @output,$_->id}}
	return \@output;
}

sub get_path_between_nodes {
	# by default, node1 is deeper than node2
	my ($node1,$node2) = @_;
#	say $node1,"\t",$node1->bootstrap,"\t",$node2,"\t",$node2->bootstrap;
	my @path = ();
	my @tmp = ();
	while (!grep {$_ eq $node2} @tmp) {
		@tmp = $node1->each_Descendent;
		foreach my $i (@tmp) {
			if (grep {$_ eq $node2} $i->get_all_Descendents) {
				push @path,$i;
				$node1 = $i;
			}
		}
	}	
#	for (@path){say $_,"\t",$_->bootstrap}
	return \@path;
}