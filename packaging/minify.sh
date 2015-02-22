#! /bin/sh

# Remove tests
rm -rf lib/app/ruby/*/gems/*/test
rm -rf lib/app/ruby/*/gems/*/tests
rm -rf lib/app/ruby/*/gems/*/spec
rm -rf lib/app/ruby/*/gems/*/features
rm -rf lib/app/ruby/*/gems/*/benchmark

# Remove documentation
rm -f lib/app/ruby/*/gems/*/README*
rm -f lib/app/ruby/*/gems/*/CHANGE*
rm -f lib/app/ruby/*/gems/*/Change*
rm -f lib/app/ruby/*/gems/*/COPYING*
rm -f lib/app/ruby/*/gems/*/LICENSE*
rm -f lib/app/ruby/*/gems/*/MIT-LICENSE*
rm -f lib/app/ruby/*/gems/*/*.txt
rm -f lib/app/ruby/*/gems/*/*.md
rm -f lib/app/ruby/*/gems/*/*.rdoc
rm -rf lib/app/ruby/*/gems/*/doc
rm -rf lib/app/ruby/*/gems/*/docs
rm -rf lib/app/ruby/*/gems/*/example
rm -rf lib/app/ruby/*/gems/*/examples
rm -rf lib/app/ruby/*/gems/*/sample
rm -rf lib/app/ruby/*/gems/*/doc-api
find lib/app/ruby -name '*.md' | xargs rm -f

# Remove misc unnecessary files
rm -rf lib/app/ruby/*/gems/*/.gitignore
rm -rf lib/app/ruby/*/gems/*/.travis.yml

# Remove leftover native extension sources
rm -f lib/app/ruby/*/gems/*/ext/Makefile
rm -f lib/app/ruby/*/gems/*/ext/*/Makefile
find lib/app/ruby -name '*.c' | xargs rm -f
find lib/app/ruby -name '*.cpp' | xargs rm -f
find lib/app/ruby -name '*.h' | xargs rm -f
find lib/app/ruby -name '*.rl' | xargs rm -f
find lib/app/ruby -name 'extconf.rb' | xargs rm -f

# Remove Java files. They're only used for JRuby support
find lib/app/ruby -name '*.java' | xargs rm -f
find lib/app/ruby -name '*.class' | xargs rm -f
