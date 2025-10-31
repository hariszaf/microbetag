import org.yaml.snakeyaml.Yaml
import groovy.yaml.YamlSlurper

def READPARAMSFILE(defaultFile) {

    // 1. Determine which params file to load
    def paramsFile = params.get('paramsFile', "${defaultFile}")

    // 2. Parse YAML file
    def yaml = new YamlSlurper().parse(new File(paramsFile))

    // 3. Assign YAML values to params (only if not already specified via CLI)
    yaml.each { k, v ->
        if (!params.containsKey(k)) {
            params[k] = v
        }
    }
    return params
}


def SANITIZE(name) {
    return name.bytes.encodeBase64().toString() //.replace('+','-').replace('/','_').replaceAll('=+$','')
}

def UNSANITIZE(name) {
    return new String(name.decodeBase64())
}

def SANITIZE_CH(pattern) {
    return Channel
        .fromPath(pattern)
        .map { file ->
            def orig_name = file.getName() 
            def safe_name = SANITIZE(orig_name)
            tuple(orig_name, safe_name, file)
        }
}

process PREP_FILES {

    tag "Copy of the input files with sanitized names"
    if (params.debug_decompress) {
        publishDir "${params.outdir}/data", mode: 'copy', overwrite: true
    }
    container "microbetag"

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

    if (params.debug_decompress) {
        publishDir "${params.outdir}/decompr", mode: 'copy', overwrite: true
    }
    container "microbetag"

    input:
    tuple val(orig_name), val(orig_name_decomp), path(file)

    output:
    // tuple val(orig_name), val(orig_name_decomp), path(orig_name_decomp)
    path(orig_name_decomp)

    script:
    """
    if file -L "${file}" | grep -q 'gzip compressed data'; then
        gunzip -c "${file}" > "${orig_name_decomp}"
    else
        cp "${file}" "${orig_name_decomp}"
    fi
    """
}


