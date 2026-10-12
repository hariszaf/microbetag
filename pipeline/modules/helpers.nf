//---------
// FUNCTIONS
//---------

def sanitize(name) {
    return name.bytes.encodeBase64().toString() //.replace('+','-').replace('/','_').replaceAll('=+$','')
}

def sanitizeChannel(fileList) {    
    return channel.fromList(fileList)
        .map { file ->
            def orig_name = file.getName() 
            def safe_name = sanitize(orig_name)
            tuple(orig_name, safe_name, file)
        }
}

def getInputFiles(String inputDir) {
    // Detects and sanitizes sequencing (FASTA) files 

    def input_files = file(inputDir)
    def faa_files   = input_files.listFiles().findAll { f -> f.name =~ /\.(faa|faa\.gz)$/ }
    def nucl_files  = input_files.listFiles().findAll { f -> f.name =~ /\.(fa|fna|fasta|fa\.gz|fna\.gz|fasta\.gz)$/ }
    
    if (faa_files && !nucl_files) {
        log.info "Detected protein FASTA files (*.faa or *.faa.gz)"
        def pattern = "${inputDir}/*.{faa,faa.gz}"
        return [
            pattern: pattern,
            is_faa: true,
            files_ch: sanitizeChannel(faa_files)
        ]
    } else if (nucl_files && !faa_files) {
        log.info "Detected nucleotide FASTA files (*.fa, *.fna, *.fasta, etc.)"
        def pattern = "${inputDir}/*.{fa,fasta,fna,fa.gz,fasta.gz,fna.gz}"
        return [
            pattern: pattern,
            is_faa: false,
            files_ch: sanitizeChannel(nucl_files)
        ]
    } else if (faa_files && nucl_files) {
        error "Mixed FASTA file types detected (both nucleotide and protein). Please separate them."
    } else {
        error "No input FASTA files found in ${inputDir}"
    }
}


def chunkFiles(fileChannel, maxForks) {
    // Split the files in at most maxForks chunks, to be processed in parallel
    return fileChannel.collect().flatMap { files ->
        log.info "Number of input files found: ${files.size()}"
        def chunk_size = Math.ceil(files.size() / maxForks) as int
        log.info "Chunk size = ${chunk_size}"
        return files.collate(chunk_size)
    }
}


//---------
// PROCESSES 
//---------

process SAFENAME_FILES {

    tag "Copy of the input files with sanitized names"
    publishDir "${params.outdir}/data", mode: 'copy', overwrite: true, enabled: params.debug_decompress ?: false
    container params.microbetag_image

    input:
    tuple val(orig_name), val(safe_name), path(file)

    output:
    tuple val(orig_name), path("${safe_name}")

    script:
    """
    # Always operate on the real file ($file)
    cp $file ${safe_name}
    chmod 644 ${safe_name}
    """
}


process GUNZIP {

    tag "Gunzip files"
    publishDir "${params.outdir}/decompr", mode: 'copy', overwrite: true, enabled: params.debug_decompress ?: false
    container params.microbetag_image

    input:
    tuple val(orig_name), val(orig_name_decomp), path(file)

    output:
    path orig_name_decomp

    script:
    """
    if file -L "${file}" | grep -q 'gzip compressed data'; then
        gunzip -c "${file}" > "${orig_name_decomp}"
    else
        cp "${file}" "${orig_name_decomp}"
    fi
    """
}

