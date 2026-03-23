#include "git.h"

// clang-format off
bool GitMetadata::Populated() { return false; }
bool GitMetadata::AnyUncommittedChanges() { return false; }
std::string GitMetadata::AuthorName() { return ""; }
std::string GitMetadata::AuthorEmail() { return ""; }
std::string GitMetadata::CommitSHA1() { return ""; }
std::string GitMetadata::CommitDate() { return ""; }
std::string GitMetadata::CommitSubject() { return ""; }
std::string GitMetadata::CommitBody() { return ""; }
std::string GitMetadata::Describe() { return "unknown"; }
// clang-format on
