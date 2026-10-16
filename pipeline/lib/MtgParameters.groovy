/*
 * The parameters of the pipeline, from parameters.yaml: the one place where they are described.
 * main.nf takes the defaults from it, checks the parameters of a run, and writes the templates and the help.
 * (Nextflow loads the classes of lib/ itself.)
 */

import groovy.yaml.YamlSlurper

class MtgParameters {

    static final int WIDTH = 110

    // ---------------- Commands, as the user runs them: through microbetag-nf (which sets MICROBETAG_NF), or not ----

    static boolean viaWrapper() { System.getenv('MICROBETAG_NF') as boolean }

    static String runCommand(String file) {
        viaWrapper() ? "microbetag-nf run ${file}" : "nextflow run main.nf -params-file ${file}"
    }

    static String helpCommand(String wf = null) {
        viaWrapper() ? "microbetag-nf describe${wf ? ' ' + wf : ''}" : "nextflow run main.nf --help${wf ? ' --workflow ' + wf : ''}"
    }

    static String writeCommand(String wf = '<workflow>', String file = '<file>.yaml') {
        viaWrapper() ? "microbetag-nf write-params ${wf} ${file}" : "nextflow run main.nf --workflow ${wf} --write-params ${file}"
    }

    final Map spec

    MtgParameters(String path) {
        spec = new YamlSlurper().parse(new File(path)) as Map
    }

    // {name: entry} of the parameters of the templates, and of the internal ones
    Map all() {
        def params = [:]
        spec.sections.each { section -> params.putAll(section.parameters as Map) }
        params.putAll((spec.internal ?: [:]) as Map)
        return params
    }

    // Those of a workflow, by section: [[title, {name: entry}], ...]
    List sectionsOf(String wf) {
        spec.sections.collect { section ->
            [section.title, (section.parameters as Map).findAll { name, p -> wf == null || wf in (p.workflows ?: []) }]
        }.findAll { it[1] }
    }

    static boolean isGroup(Map p) { p.containsKey('parameters') }

    // The default of a parameter; of a group, the defaults of its parameters
    static Object defaultOf(Map p) {
        isGroup(p) ? (p.parameters as Map).collectEntries { k, v -> [k, defaultOf(v as Map)] } : p.value
    }

    static Object templateOf(Map p) {
        isGroup(p) ? (p.parameters as Map).collectEntries { k, v -> [k, templateOf(v as Map)] }
                   : (p.containsKey('template') ? p.template : p.value)
    }

    // The defaults go to the parameters the run leaves out (nextflow.config already sets them; this also covers a run
    // without it); a group given in part keeps the defaults of the rest
    void applyDefaults(Map params) {
        all().each { name, p ->
            def d = defaultOf(p as Map)
            if (!params.containsKey(name)) {
                params[name] = d
            } else if (isGroup(p as Map)) {
                params[name] = merge(d as Map, params[name] instanceof Map ? params[name] as Map : [:])
            }
        }
        if (params.recon_files == null) {
            params.recon_files = params.faa_dir ?: params.genomes
        }
    }

    static Map merge(Map defaults, Map given) {
        def out = new LinkedHashMap(defaults)
        given.each { k, v ->
            out[k] = (v instanceof Map && defaults[k] instanceof Map) ? merge(defaults[k] as Map, v as Map) : v
        }
        return out
    }

    // Parameters the pipeline does not know (e.g. a typo), also within groups: ["flashweave_args.max_K", ...]
    List unknown(Map params) {
        def found = []
        def known = all()
        params.each { name, value ->
            if (!known.containsKey(name)) {
                found << name
            } else if (isGroup(known[name] as Map) && value instanceof Map) {
                found.addAll(unknownIn(known[name] as Map, value as Map, name as String))
            }
        }
        return found.sort()
    }

    static List unknownIn(Map group, Map value, String prefix) {
        def found = []
        value.each { k, v ->
            def p = (group.parameters as Map)[k]
            if (p == null) {
                found << "${prefix}.${k}".toString()
            } else if (isGroup(p as Map) && v instanceof Map) {
                found.addAll(unknownIn(p as Map, v as Map, "${prefix}.${k}".toString()))
            }
        }
        return found
    }

    // Required parameters of a workflow that the run does not set
    List missing(Map params, String wf) {
        all().findAll { name, p ->
            def req = p.required
            (req == true || (req instanceof List && wf in req)) && (params[name] == null || params[name] == '')
        }.keySet() as List
    }

    // ---------------- Templates ----------------

    // The template of a workflow, to be written to `file`; it sets the workflow, so that the file runs on its own
    String template(String wf, String file = "${wf}.yaml") {
        def rule = '# ' + '-' * (WIDTH - 2)
        def fileName = new File(file).name
        def intro = (spec.intro as String).replace('{workflow}', wf).replace('{file}', fileName)
            .replace('{run}', "microbetag-nf run ${fileName}  (or: nextflow run main.nf -params-file ${fileName})")
            .replace('{help}', "microbetag-nf describe [PARAMETER] ${wf}  (or: nextflow run main.nf --help [PARAMETER] --workflow ${wf})")
        def out = ["# ${spec.title}, ${wf} workflow: ${spec.workflows[wf]}".toString(), '#', comment(intro), '',
                   '# The workflow this file runs (precalc, annotate, prep or network).', "workflow: ${wf}".toString(), '']
        sectionsOf(wf).each { title, params ->
            out += [rule, "# ${title}".toString(), rule, '']
            params.each { name, p -> out += entry(name as String, p as Map, 0) + [''] }
        }
        return out.join('\n')
    }

    static List entry(String name, Map p, int indent) {
        def pad  = ' ' * indent
        def desc = p.description as String
        if (p.required) desc += ' Required.'
        if (p.example != null) desc += " E.g. ${p.example}"
        def lines = comment(desc, WIDTH - indent).readLines().collect { pad + it }
        if (isGroup(p)) {
            lines << "${pad}${name}:".toString()
            (p.parameters as Map).each { k, v -> lines.addAll(entry(k as String, v as Map, indent + 2)) }
        } else {
            lines << "${pad}${name}: ${scalar(templateOf(p))}".toString().replaceAll(/\s+$/, '')
        }
        return lines
    }

    static String comment(String text, int width = WIDTH) {
        def lines = []
        text.trim().split(/\n\s*\n/).eachWithIndex { paragraph, i ->
            if (i > 0) lines << ''
            lines.addAll(wrap(paragraph.split(/\s+/).join(' '), width - 2))
        }
        return lines.collect { "# ${it}".toString().replaceAll(/\s+$/, '') }.join('\n')
    }

    static List wrap(String text, int width) {
        def lines = []
        def line = ''
        text.split(' ').each { word ->
            if (line && (line.size() + 1 + word.size()) > width) {
                lines << line
                line = word
            } else {
                line = line ? "${line} ${word}" : word
            }
        }
        if (line) lines << line
        return lines
    }

    // A YAML scalar: plain when it reads back as the same string, else quoted
    static String scalar(Object v) {
        if (v == null) return ''
        if (v instanceof Boolean || v instanceof Number) return v.toString()
        def s = v.toString()
        def plain = s && !(s ==~ /(?i)(true|false|yes|no|on|off|null|~|[-+]?[0-9][0-9._eE+-]*)/) &&
                    !(s =~ /[:#\[\]{},&*!|>'"%@`\t\n]|^[-? ]|\s$/)
        return plain ? s : groovy.json.JsonOutput.toJson(s)
    }

    // The parameters of a run as YAML (e.g. for microbetag-annotate, which reads them)
    static String toYaml(Map values, int indent = 0) {
        def pad = ' ' * indent
        values.collect { k, v ->
            v instanceof Map && v ? "${pad}${k}:\n" + toYaml(v as Map, indent + 2)
                                  : "${pad}${k}: ${v instanceof Map ? '{}' : scalar(v)}".replaceAll(/\s+$/, '')
        }.join('\n')
    }

    // ---------------- Help ----------------

    String describe(String name, String wf) {
        def params = all()
        if (name) {
            def path = name.tokenize('.')
            def p = params[path[0]] as Map
            path.drop(1).each { k -> p = p && isGroup(p) ? (p.parameters as Map)[k] as Map : null }
            if (p == null) {
                throw new IllegalArgumentException("Unknown parameter: ${name}. See: ${helpCommand()}")
            }
            return describeOne(name, p, '')
        }
        def out = []
        if (!wf) {
            def describeOne = viaWrapper() ? 'microbetag-nf describe <parameter>' : 'nextflow run main.nf --help <parameter>'
            out += ['microbetag pipeline: from your genomes to a microbetag-annotated co-occurrence network', '',
                    'Workflows:']
            spec.workflows.each { k, v -> out << "  ${k.padRight(10)} ${v}".toString() }
            out += ['', 'A run, step by step:',
                    "  1. write the params file of a workflow:  ${writeCommand('precalc', 'my_precalc.yaml')}".toString(),
                    '  2. fill it in: your data, the steps you want; a parameter left out takes its default',
                    "  3. run it:                                ${runCommand('my_precalc.yaml')}".toString(),
                    '',
                    "Describe a parameter:                     ${describeOne}".toString(),
                    "Only the parameters of a workflow:        ${helpCommand('<workflow>')}".toString(),
                    '', '']
        }
        out += [wf ? "Parameters of the ${wf} workflow (${spec.workflows[wf]}):" : 'All the parameters:', '']
        sectionsOf(wf).each { title, ps ->
            out << "${title}".toString()
            ps.each { n, p -> out << describeOne(n as String, p as Map, '  ') }
            out << ''
        }
        out << "Write the params file of a workflow, to fill in: ${writeCommand()}".toString()
        return out.join('\n')
    }

    static String describeOne(String name, Map p, String pad) {
        def info = isGroup(p) ? 'group' : "default: ${p.value == null ? 'not set' : p.value}"
        if (p.workflows) info += "; ${(p.workflows as List).join(', ')}"
        if (p.required) info += '; required'
        def lines = ["${pad}${name}  (${info})".toString()]
        lines.addAll(wrap((p.description as String).split(/\s+/).join(' '), WIDTH - pad.size() - 4).collect { "${pad}    ${it}".toString() })
        if (isGroup(p)) {
            (p.parameters as Map).each { k, v -> lines << describeOne("${name}.${k}".toString(), v as Map, pad + '  ') }
        }
        return lines.join('\n')
    }
}
