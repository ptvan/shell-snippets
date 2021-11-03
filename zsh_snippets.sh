##### SET SHELL OPTIONS
setopt 

# mmv file*.dat file*.txt
autoload -U zmv
alias mmv="noglob zmv -W"

# `take` will recursively mkdir and cd to deepest level
take /some/really/deeply/nested/directory
