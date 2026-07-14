# Security policy

MPPH is a command-line research tool. It makes read-only HTTPS requests to the
public KEGG REST API and writes files (figures, tables, a cache) under the
directories you specify. It does not require credentials and does not transmit
your data anywhere except to KEGG when you query it.

## Reporting a vulnerability

If you find a security issue (for example, a way this tool could be made to
write outside its output directory, or a problem in how cached responses are
handled), please open a GitHub issue marked "security", or contact the
maintainer privately via their GitHub profile. Please do not include sensitive
details in a public issue.
