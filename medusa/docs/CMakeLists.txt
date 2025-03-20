# generate docs
add_custom_target(docs COMMAND doxygen docs/Doxyfile > /dev/null
                  WORKING_DIRECTORY ${MEDUSA_ROOT})
# deploy docs
add_custom_target(deploy COMMAND ./scripts/deploy.sh DEPENDS docs
                  WORKING_DIRECTORY ${MEDUSA_ROOT})
# backup wiki
add_custom_target(static_wiki COMMAND ./scripts/backup_wiki_static.sh
                  WORKING_DIRECTORY ${MEDUSA_ROOT})
