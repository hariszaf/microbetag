import org.yaml.snakeyaml.Yaml
import groovy.yaml.YamlSlurper

def readParamsFile(defaultFile) {

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
