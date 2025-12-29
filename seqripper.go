// SeqRipper - Extract genes from GenBank files
// A pure Go tool for extracting specific genes from local GenBank (.gb) files

package main

import (
	"bufio"
	"flag"
	"fmt"
	"os"
	"path/filepath"
	"regexp"
	"sort"
	"strconv"
	"strings"
	"sync"
)

const version = "2.0.0"

// Compiled regex patterns for location parsing (compiled once at startup)
var (
	rangePattern  = regexp.MustCompile(`<?(\d+)\.\.>?(\d+)`)
	singlePattern = regexp.MustCompile(`^<?(\d+)>?$`)
)

// Gene synonyms for matching annotations
var geneSynonyms = map[string][]string{
	"12S":   {"12S", "12S rRNA", "s-rRNA", "12S ribosomal RNA", "rrnS"},
	"16S":   {"16S", "16S rRNA", "l-rRNA", "16S ribosomal RNA", "rrnL"},
	"COI":   {"COI", "COXI", "COX1", "CO1", "cytochrome oxidase subunit 1", "cytochrome c oxidase subunit I", "cytochrome oxidase subunit I"},
	"COII":  {"COII", "COXII", "COX2", "CO2", "cytochrome oxidase subunit 2", "cytochrome c oxidase subunit II"},
	"COIII": {"COIII", "COXIII", "COX3", "CO3", "cytochrome oxidase subunit 3"},
	"CYTB":  {"CYTB", "Cyt b", "cytochrome b", "cob", "CYB"},
	"ND1":   {"ND1", "NAD1", "NADH1", "NADH dehydrogenase subunit 1"},
	"ND2":   {"ND2", "NAD2", "NADH2", "NADH dehydrogenase subunit 2"},
	"ND3":   {"ND3", "NAD3", "NADH3", "NADH dehydrogenase subunit 3"},
	"ND4":   {"ND4", "NAD4", "NADH4", "NADH dehydrogenase subunit 4"},
	"ND4L":  {"ND4L", "NAD4L", "NADH4L", "NADH dehydrogenase subunit 4L"},
	"ND5":   {"ND5", "NAD5", "NADH5", "NADH dehydrogenase subunit 5"},
	"ND6":   {"ND6", "NAD6", "NADH6", "NADH dehydrogenase subunit 6"},
	"ATP6":  {"ATP6", "ATPase6", "ATPase 6", "ATP synthase F0 subunit 6"},
	"ATP8":  {"ATP8", "ATPase8", "ATPase 8", "ATP synthase F0 subunit 8"},
}

var validFeatureTypes = map[string]bool{
	"rRNA":    true,
	"CDS":     true,
	"gene":    true,
	"tRNA":    true,
	"misc_RNA": true,
}

// GenBank record structures
type Location struct {
	Start      int
	End        int
	Complement bool
}

type Feature struct {
	Type       string
	Locations  []Location
	Qualifiers map[string]string
	Complement bool
}

type GenBankRecord struct {
	Accession  string
	Organism   string
	Sequence   string
	Features   []Feature
}

// Extraction result for logging
type ExtractionResult struct {
	Accession string
	Gene      string
	LengthBP  int
	Status    string
	Filename  string
}

// Command line arguments
type Config struct {
	InputFiles []string
	Genes      []string
	OutputDir  string
	LogFile    string
	ExtractAll bool
	SkipTRNA   bool
	Quiet      bool
	Verbose    bool
}

// Custom string slice flag type for comma-separated values
type stringSlice []string

func (s *stringSlice) String() string {
	return strings.Join(*s, ", ")
}

func (s *stringSlice) Set(value string) error {
	parts := strings.Split(value, ",")
	for _, p := range parts {
		p = strings.TrimSpace(p)
		if p != "" {
			*s = append(*s, p)
		}
	}
	return nil
}

// Track used filenames for duplicate handling
type FilenameTracker struct {
	mu    sync.Mutex
	used  map[string]int
}

func NewFilenameTracker() *FilenameTracker {
	return &FilenameTracker{
		used: make(map[string]int),
	}
}

func (ft *FilenameTracker) GetUniqueFilename(base string) (string, bool) {
	ft.mu.Lock()
	defer ft.mu.Unlock()
	
	count, exists := ft.used[base]
	if !exists {
		ft.used[base] = 1
		return base, false
	}
	
	// Duplicate - create numbered version
	ft.used[base] = count + 1
	ext := filepath.Ext(base)
	name := strings.TrimSuffix(base, ext)
	return fmt.Sprintf("%s_%d%s", name, count, ext), true
}

// Parse GenBank file
func parseGenBank(filename string) (*GenBankRecord, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, fmt.Errorf("failed to open file: %w", err)
	}
	defer file.Close()

	record := &GenBankRecord{
		Features: make([]Feature, 0),
	}

	scanner := bufio.NewScanner(file)
	// Increase buffer size for long lines
	buf := make([]byte, 0, 64*1024)
	scanner.Buffer(buf, 1024*1024)

	var inFeatures, inOrigin bool
	var currentFeature *Feature
	var featureLines []string
	var seqBuilder strings.Builder

	for scanner.Scan() {
		line := scanner.Text()

		// VERSION line - get accession
		if strings.HasPrefix(line, "VERSION") {
			fields := strings.Fields(line)
			if len(fields) >= 2 {
				record.Accession = fields[1]
			}
			continue
		}

		// ORGANISM in source section
		if strings.Contains(line, "ORGANISM") && record.Organism == "" {
			// ORGANISM is typically followed by the name on the same line
			parts := strings.SplitN(line, "ORGANISM", 2)
			if len(parts) == 2 {
				record.Organism = strings.TrimSpace(parts[1])
			}
			continue
		}

		// Section markers
		if strings.HasPrefix(line, "FEATURES") {
			inFeatures = true
			continue
		}
		if strings.HasPrefix(line, "ORIGIN") {
			inFeatures = false
			inOrigin = true
			// Save any pending feature
			if currentFeature != nil && len(featureLines) > 0 {
				parseFeatureBlock(currentFeature, featureLines)
				record.Features = append(record.Features, *currentFeature)
			}
			continue
		}
		if strings.HasPrefix(line, "//") {
			break
		}

		// Parse features section
		if inFeatures {
			if len(line) > 5 && line[5] != ' ' {
				// New feature type at position 5
				// Save previous feature
				if currentFeature != nil && len(featureLines) > 0 {
					parseFeatureBlock(currentFeature, featureLines)
					record.Features = append(record.Features, *currentFeature)
				}
				
				// Start new feature
				fields := strings.Fields(line)
				if len(fields) >= 2 {
					currentFeature = &Feature{
						Type:       fields[0],
						Qualifiers: make(map[string]string),
						Locations:  make([]Location, 0),
					}
					featureLines = []string{line}
				}
			} else if currentFeature != nil {
				// Continuation of current feature
				featureLines = append(featureLines, line)
			}
		}

		// Parse sequence in ORIGIN section
		if inOrigin {
			// Remove numbers and spaces, keep only sequence
			for _, ch := range line {
				if ch >= 'a' && ch <= 'z' {
					seqBuilder.WriteByte(byte(ch - 32)) // Convert to uppercase
				} else if ch >= 'A' && ch <= 'Z' {
					seqBuilder.WriteByte(byte(ch))
				}
			}
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading file: %w", err)
	}

	record.Sequence = seqBuilder.String()

	// Fallback for accession if VERSION not found
	if record.Accession == "" {
		record.Accession = strings.TrimSuffix(filepath.Base(filename), filepath.Ext(filename))
	}

	return record, nil
}

// Parse a feature block (type line + qualifiers)
func parseFeatureBlock(feature *Feature, lines []string) {
	if len(lines) == 0 {
		return
	}

	// First line has the location
	fields := strings.Fields(lines[0])
	if len(fields) >= 2 {
		parseLocation(feature, fields[1])
	}

	// Remaining lines have qualifiers
	var currentQualifier string
	var currentValue strings.Builder

	for i := 1; i < len(lines); i++ {
		line := strings.TrimSpace(lines[i])
		
		if strings.HasPrefix(line, "/") {
			// Save previous qualifier
			if currentQualifier != "" {
				feature.Qualifiers[currentQualifier] = strings.Trim(currentValue.String(), "\"")
			}
			
			// Parse new qualifier
			line = strings.TrimPrefix(line, "/")
			if idx := strings.Index(line, "="); idx != -1 {
				currentQualifier = line[:idx]
				currentValue.Reset()
				currentValue.WriteString(strings.TrimPrefix(line[idx+1:], "\""))
			} else {
				currentQualifier = line
				currentValue.Reset()
			}
		} else if currentQualifier != "" {
			// Continuation of qualifier value
			currentValue.WriteString(strings.Trim(line, "\""))
		}
	}
	
	// Save last qualifier
	if currentQualifier != "" {
		feature.Qualifiers[currentQualifier] = strings.Trim(currentValue.String(), "\"")
	}
}

// Parse location string (handles complement, join, simple ranges)
func parseLocation(feature *Feature, locStr string) {
	// Check for complement
	if strings.HasPrefix(locStr, "complement(") {
		feature.Complement = true
		locStr = strings.TrimPrefix(locStr, "complement(")
		locStr = strings.TrimSuffix(locStr, ")")
	}

	// Check for join
	if strings.HasPrefix(locStr, "join(") {
		locStr = strings.TrimPrefix(locStr, "join(")
		locStr = strings.TrimSuffix(locStr, ")")
	}

	// Parse individual ranges
	matches := rangePattern.FindAllStringSubmatch(locStr, -1)

	for _, match := range matches {
		if len(match) >= 3 {
			start, _ := strconv.Atoi(match[1])
			end, _ := strconv.Atoi(match[2])
			feature.Locations = append(feature.Locations, Location{
				Start: start,
				End:   end,
			})
		}
	}

	// Handle single position (rare but possible)
	if len(feature.Locations) == 0 {
		if match := singlePattern.FindStringSubmatch(locStr); len(match) >= 2 {
			pos, _ := strconv.Atoi(match[1])
			feature.Locations = append(feature.Locations, Location{
				Start: pos,
				End:   pos,
			})
		}
	}
}

// Extract sequence for a feature
func extractSequence(record *GenBankRecord, feature *Feature) string {
	var seq strings.Builder
	
	for _, loc := range feature.Locations {
		start := loc.Start - 1 // Convert to 0-indexed
		end := loc.End
		
		if start < 0 {
			start = 0
		}
		if end > len(record.Sequence) {
			end = len(record.Sequence)
		}
		if start >= len(record.Sequence) {
			continue
		}
		
		seq.WriteString(record.Sequence[start:end])
	}

	result := seq.String()
	
	if feature.Complement {
		result = reverseComplement(result)
	}
	
	return result
}

// Reverse complement a DNA sequence
func reverseComplement(seq string) string {
	complement := map[byte]byte{
		'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
		'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
		'N': 'N', 'n': 'n',
		'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
		'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
		'D': 'H', 'H': 'D',
	}
	
	result := make([]byte, len(seq))
	for i := 0; i < len(seq); i++ {
		c := seq[len(seq)-1-i]
		if comp, ok := complement[c]; ok {
			result[i] = comp
		} else {
			result[i] = c
		}
	}
	return string(result)
}

// Get gene name from feature
func getGeneName(feature *Feature) string {
	if gene, ok := feature.Qualifiers["gene"]; ok {
		return gene
	}
	if product, ok := feature.Qualifiers["product"]; ok {
		return product
	}
	return ""
}

// Check if feature is a tRNA (by type or name)
func isTRNA(feature *Feature, geneName string) bool {
	if feature.Type == "tRNA" {
		return true
	}
	geneNameLower := strings.ToLower(geneName)
	return strings.HasPrefix(geneNameLower, "trn") || strings.Contains(geneNameLower, "trna")
}

// Match annotation to target gene using synonyms
func matchGene(annotation string, targetGenes map[string][]string) string {
	annotationLower := strings.ToLower(strings.TrimSpace(annotation))
	
	for gene, synonyms := range targetGenes {
		// Exact match first
		for _, syn := range synonyms {
			if annotationLower == strings.ToLower(syn) {
				return gene
			}
		}
		// Partial match
		for _, syn := range synonyms {
			if strings.Contains(annotationLower, strings.ToLower(syn)) {
				return gene
			}
		}
	}
	return ""
}

// Build target genes map from user input
func parseTargetGenes(genes []string) map[string][]string {
	result := make(map[string][]string)
	
	for _, gene := range genes {
		geneUpper := strings.ToUpper(gene)
		if synonyms, ok := geneSynonyms[geneUpper]; ok {
			result[geneUpper] = synonyms
		} else {
			// Custom gene - basic synonyms
			result[geneUpper] = []string{gene, strings.ToUpper(gene), strings.ToLower(gene)}
		}
	}
	return result
}

// Clean name for use in filenames
func cleanName(name string) string {
	name = strings.ReplaceAll(name, " ", "_")
	name = strings.ReplaceAll(name, "/", "_")
	name = strings.ReplaceAll(name, "\\", "_")
	name = strings.ReplaceAll(name, ":", "_")
	return name
}

// Write FASTA file
func writeFASTA(filename, header, sequence string) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	writer := bufio.NewWriter(file)
	fmt.Fprintf(writer, ">%s\n", header)
	
	// Write sequence in 60-character lines
	for i := 0; i < len(sequence); i += 60 {
		end := i + 60
		if end > len(sequence) {
			end = len(sequence)
		}
		fmt.Fprintf(writer, "%s\n", sequence[i:end])
	}
	
	return writer.Flush()
}

// Get location string for FASTA header
func getLocationString(feature *Feature) string {
	if len(feature.Locations) == 0 {
		return "1-1"
	}
	
	// Get overall start and end
	start := feature.Locations[0].Start
	end := feature.Locations[len(feature.Locations)-1].End
	
	return fmt.Sprintf("%d-%d", start, end)
}

// Process a single GenBank file
func processFile(filename string, config *Config, targetGenes map[string][]string, 
	tracker *FilenameTracker, resultsChan chan<- ExtractionResult, verbose bool) {
	
	record, err := parseGenBank(filename)
	if err != nil {
		if verbose {
			fmt.Fprintf(os.Stderr, "  ✗ Error parsing %s: %v\n", filename, err)
		}
		resultsChan <- ExtractionResult{
			Accession: filepath.Base(filename),
			Gene:      "-",
			LengthBP:  0,
			Status:    "error",
			Filename:  "-",
		}
		return
	}

	if verbose {
		fmt.Fprintf(os.Stderr, "  Processing: %s (%s)\n", record.Accession, record.Organism)
	}

	foundGenes := make(map[string]bool)

	for _, feature := range record.Features {
		// Check feature type
		if !validFeatureTypes[feature.Type] {
			continue
		}

		geneName := getGeneName(&feature)
		if geneName == "" {
			continue
		}

		// Skip tRNA if requested
		if config.SkipTRNA && isTRNA(&feature, geneName) {
			continue
		}

		var targetName string
		if config.ExtractAll {
			targetName = cleanName(geneName)
		} else {
			targetName = matchGene(geneName, targetGenes)
		}

		if targetName == "" || foundGenes[targetName] {
			continue
		}

		// Extract sequence
		sequence := extractSequence(record, &feature)
		if len(sequence) == 0 {
			continue
		}

		// Generate filename
		baseFilename := fmt.Sprintf("%s_%s.fasta", record.Accession, targetName)
		uniqueFilename, isDuplicate := tracker.GetUniqueFilename(baseFilename)
		
		if isDuplicate && verbose {
			fmt.Fprintf(os.Stderr, "    ⚠ Duplicate: %s → %s\n", baseFilename, uniqueFilename)
		}

		outputPath := filepath.Join(config.OutputDir, uniqueFilename)

		// Create FASTA header
		locStr := getLocationString(&feature)
		header := fmt.Sprintf("%s:%s %s [gene=%s]", record.Accession, locStr, record.Organism, targetName)

		// Write file
		if err := writeFASTA(outputPath, header, sequence); err != nil {
			if verbose {
				fmt.Fprintf(os.Stderr, "    ✗ Failed to write %s: %v\n", uniqueFilename, err)
			}
			resultsChan <- ExtractionResult{
				Accession: record.Accession,
				Gene:      targetName,
				LengthBP:  0,
				Status:    "error",
				Filename:  "-",
			}
			continue
		}

		foundGenes[targetName] = true
		
		resultsChan <- ExtractionResult{
			Accession: record.Accession,
			Gene:      targetName,
			LengthBP:  len(sequence),
			Status:    "extracted",
			Filename:  uniqueFilename,
		}

		if verbose {
			fmt.Fprintf(os.Stderr, "    ✓ %s → %s (%d bp)\n", targetName, uniqueFilename, len(sequence))
		}
	}

	// Report missing genes (only in specific gene mode)
	if !config.ExtractAll {
		for gene := range targetGenes {
			if !foundGenes[gene] {
				resultsChan <- ExtractionResult{
					Accession: record.Accession,
					Gene:      gene,
					LengthBP:  0,
					Status:    "missing",
					Filename:  "-",
				}
				if verbose {
					fmt.Fprintf(os.Stderr, "    ✗ %s: not found\n", gene)
				}
			}
		}
	}
}

// Write TSV log file
func writeLog(filename string, results []ExtractionResult) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	writer := bufio.NewWriter(file)
	
	// Header
	fmt.Fprintln(writer, "Accession\tGene\tLength_bp\tStatus\tFilename")
	
	// Sort results for consistent output
	sort.Slice(results, func(i, j int) bool {
		if results[i].Accession != results[j].Accession {
			return results[i].Accession < results[j].Accession
		}
		return results[i].Gene < results[j].Gene
	})

	for _, r := range results {
		fmt.Fprintf(writer, "%s\t%s\t%d\t%s\t%s\n", 
			r.Accession, r.Gene, r.LengthBP, r.Status, r.Filename)
	}

	return writer.Flush()
}

func printUsage() {
	fmt.Fprintf(os.Stderr, `SeqRipper v%s - Extract genes from GenBank files

Usage:
  seqripper [options] <input.gb ...> <output_dir>

Gene Selection (one required):
  -g GENES      Comma-separated genes to extract (e.g., COI,12S,CYTB)
  --all         Extract ALL annotated genes

Options:
  --no-trna     Skip tRNA genes (only with --all)
  -v, --verbose Show detailed per-file output
  -q, --quiet   Suppress all output
  -h, --help    Show this help
  --version     Print version and exit

Arguments:
  input.gb      One or more GenBank files (shell globs work)
  output_dir    Output directory (created if needed)

Output:
  <output_dir>/<accession>_<gene>.fasta
  <output_dir>/seqripper_log.tsv

Examples:
  seqripper -g COI,12S,CYTB mitogenomes/*.gb results/
  seqripper --all data/*.gb extracted_genes/
  seqripper --all --no-trna -q *.gb output/

Supported genes with synonyms:
  12S, 16S, COI, COII, COIII, CYTB, ND1-ND6, ND4L, ATP6, ATP8
`, version)
}

func main() {
	// Custom flag set for better control
	flags := flag.NewFlagSet("seqripper", flag.ExitOnError)
	flags.Usage = printUsage

	// Define flags
	var genes stringSlice
	flags.Var(&genes, "g", "Genes to extract")
	
	extractAll := flags.Bool("all", false, "Extract all genes")
	skipTRNA := flags.Bool("no-trna", false, "Skip tRNA genes")
	quiet := flags.Bool("q", false, "Quiet mode")
	quietLong := flags.Bool("quiet", false, "Quiet mode")
	verbose := flags.Bool("v", false, "Verbose output")
	verboseLong := flags.Bool("verbose", false, "Verbose output")
	showVersion := flags.Bool("version", false, "Show version")
	help := flags.Bool("h", false, "Show help")
	helpLong := flags.Bool("help", false, "Show help")

	// Parse flags
	if err := flags.Parse(os.Args[1:]); err != nil {
		os.Exit(1)
	}

	// Handle help and version first
	if *help || *helpLong {
		printUsage()
		os.Exit(0)
	}

	if *showVersion {
		fmt.Printf("SeqRipper v%s\n", version)
		os.Exit(0)
	}

	// Combine quiet and verbose flags
	*quiet = *quiet || *quietLong
	*verbose = *verbose || *verboseLong

	// Positional arguments: INPUT... OUTPUT_DIR
	args := flags.Args()
	if len(args) < 2 {
		fmt.Fprintln(os.Stderr, "Error: need at least one input file and an output directory")
		printUsage()
		os.Exit(1)
	}

	// Last arg is output directory, rest are input files
	outputDir := args[len(args)-1]
	inputFiles := args[:len(args)-1]
	logFile := filepath.Join(outputDir, "seqripper_log.tsv")

	// Validate gene selection
	if !*extractAll && len(genes) == 0 {
		fmt.Fprintln(os.Stderr, "Error: either -g (genes) or --all is required")
		printUsage()
		os.Exit(1)
	}

	if *skipTRNA && !*extractAll {
		fmt.Fprintln(os.Stderr, "Warning: --no-trna only works with --all mode")
	}

	// Verify input files exist
	var validFiles []string
	for _, f := range inputFiles {
		if _, err := os.Stat(f); err == nil {
			validFiles = append(validFiles, f)
		} else if !*quiet {
			fmt.Fprintf(os.Stderr, "Warning: file not found: %s\n", f)
		}
	}

	if len(validFiles) == 0 {
		fmt.Fprintln(os.Stderr, "Error: no valid input files found")
		os.Exit(1)
	}

	// Create output directory
	if err := os.MkdirAll(outputDir, 0755); err != nil {
		fmt.Fprintf(os.Stderr, "Error creating output directory: %v\n", err)
		os.Exit(1)
	}

	// Build config
	config := &Config{
		InputFiles: validFiles,
		Genes:      genes,
		OutputDir:  outputDir,
		LogFile:    logFile,
		ExtractAll: *extractAll,
		SkipTRNA:   *skipTRNA,
		Quiet:      *quiet,
		Verbose:    *verbose,
	}

	// Build target genes map
	var targetGenes map[string][]string
	if !*extractAll {
		targetGenes = parseTargetGenes(genes)
	} else {
		targetGenes = make(map[string][]string)
	}

	// Process files concurrently
	tracker := NewFilenameTracker()
	resultsChan := make(chan ExtractionResult, 100)
	var wg sync.WaitGroup

	// Collect results in background
	var allResults []ExtractionResult
	var resultsMu sync.Mutex
	done := make(chan bool)

	go func() {
		for result := range resultsChan {
			resultsMu.Lock()
			allResults = append(allResults, result)
			resultsMu.Unlock()
		}
		done <- true
	}()

	// Process files
	for _, filename := range validFiles {
		wg.Add(1)
		go func(f string) {
			defer wg.Done()
			processFile(f, config, targetGenes, tracker, resultsChan, *verbose)
		}(filename)
	}

	wg.Wait()
	close(resultsChan)
	<-done

	// Write log file
	if err := writeLog(logFile, allResults); err != nil {
		fmt.Fprintf(os.Stderr, "Error writing log file: %v\n", err)
		os.Exit(1)
	}

	// Print summary
	if !*quiet {
		extracted := 0
		errors := 0
		for _, r := range allResults {
			switch r.Status {
			case "extracted":
				extracted++
			case "error":
				errors++
			}
		}

		fmt.Fprintln(os.Stderr, "SUMMARY:")
		fmt.Fprintf(os.Stderr, "Sequences extracted: %d\n", extracted)
		if errors > 0 {
			fmt.Fprintf(os.Stderr, "Errors: %d\n", errors)
		}
		fmt.Fprintf(os.Stderr, "Output written to: %s\n", outputDir)
	}
}
