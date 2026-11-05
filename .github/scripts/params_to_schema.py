# Usage: python .github/scripts/params_to_schema.py main.nf nextflow_schema.json
"""Generate nextflow_schema.json from the typed params block in main.nf (documentation / launch forms)."""
import json, re, sys
src = open(sys.argv[1]).read()
block = src[src.index("params {") + len("params {"): src.index("\n}\n", src.index("params {"))]
groups, cur = [], None
for line in block.splitlines():
    s = line.strip()
    if not s:
        continue
    if s.startswith("//") and ":" not in s.split("//", 1)[1][:40] or (s.startswith("//") and not cur):
        title = s.lstrip("/ ").split(":")[0]
        cur = {"title": title, "props": {}, "desc": s.lstrip("/ ")}
        groups.append(cur); continue
    if s.startswith("//"):
        title = s.lstrip("/ ").split(":")[0]
        cur = {"title": title, "props": {}, "desc": s.lstrip("/ ")}
        groups.append(cur); continue
    m = re.match(r"(\w+):\s*([\w<>]+)(\?)?\s*(?:=\s*([^/]+?))?\s*(?://\s*(.*))?$", s)
    if not m:
        raise SystemExit(f"cannot parse: {s}")
    name, typ, nullable, default, comment = m.groups()
    prop = {}
    if typ == "Boolean":
        prop["type"] = "boolean"
    elif typ == "Integer":
        prop["type"] = "integer"
    elif typ == "Float":
        prop["type"] = "number"
    elif typ == "String":
        prop["type"] = "string"
    elif typ == "Path":
        prop["type"] = "string"; prop["format"] = "path"
    elif typ.startswith("Channel<"):
        prop["type"] = "string"; prop["format"] = "file-path"; prop["mimetype"] = "text/csv"; prop["pattern"] = r"^\S+\.csv$"
    if default is not None:
        d = default.strip()
        if d in ("true", "false"):
            prop["default"] = d == "true"
        elif re.fullmatch(r"-?\d+", d):
            prop["default"] = int(d)
        elif re.fullmatch(r"-?\d+\.\d*", d):
            prop["default"] = float(d)
        else:
            prop["default"] = d.strip("'\"")
    prop["description"] = (comment or name.replace("_", " ")).strip()
    cur["props"][name] = (prop, nullable is None and default is None)
defs, allof = {}, []
for g in groups:
    key = re.sub(r"[^a-z0-9]+", "_", g["title"].lower()).strip("_")
    d = {"title": g["title"], "type": "object", "description": g["desc"], "properties": {k: v[0] for k, v in g["props"].items()}}
    req = [k for k, v in g["props"].items() if v[1]]
    if req:
        d["required"] = req
    defs[key] = d; allof.append({"$ref": f"#/$defs/{key}"})
defs["output_options"] = {"title": "Output options", "type": "object", "description": "Configuration-level options (nextflow.config)",
    "properties": {"outdir": {"type": "string", "format": "directory-path", "default": "results", "description": "Output directory"},
                   "publish_dir_mode": {"type": "string", "default": "copy", "description": "How outputs are published",
                                        "enum": ["symlink", "rellink", "link", "copy", "copyNoFollow", "move"]}}}
allof.append({"$ref": "#/$defs/output_options"})
schema = {"$schema": "https://json-schema.org/draft/2020-12/schema",
          "$id": "https://raw.githubusercontent.com/barbarahelena/gutsnv/main/nextflow_schema.json",
          "title": "gutsnv pipeline parameters",
          "description": "Generated from the typed params block in main.nf. Parameter validation is done by Nextflow itself; this schema is for documentation and launch forms.",
          "type": "object", "$defs": defs, "allOf": allof}
json.dump(schema, open(sys.argv[2], "w"), indent=4); open(sys.argv[2], "a").write("\n")
print(", ".join(f"{g['title']} ({len(g['props'])})" for g in groups))
