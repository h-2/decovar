
#include <iostream>

#include <sharg/all.hpp>

#include "allele/allele.hpp"
#include "binalleles/binalleles.hpp"
#include "misc.hpp"

int main(int argc, char ** argv)
{
    std::ios::sync_with_stdio(false);

    sharg::parser parser{
      "decovar",
      argc,
      argv,
      sharg::update_notifications::off,
      {"allele", "binalleles"}
    };
    parser.info.author            = "Hannes Hauswedell";
    parser.info.short_description = "deCODE variant tools.";
    parser.info.version           = "0.1.0";
    //TODO add description

#ifdef NDEBUG
    try
    {
#endif
        parser.parse();
        sharg::parser & sub_parser = parser.get_sub_parser();

        if (sub_parser.info.app_name == std::string_view{"decovar-allele"})
            allele(sub_parser);
        else if (sub_parser.info.app_name == std::string_view{"decovar-binalleles"})
            _binalleles::main(sub_parser);
        else
            throw decovar_error{"Unhandled subcommand {} encountered. ", sub_parser.info.app_name};
#ifdef NDEBUG
    }
    catch (sharg::parser_error const & ext)
    {
        fmt::print(stderr, "[Parsing error] {}\n", ext.what());
        return -1;
    }
    catch (bio::io::bio_error const & ext)
    {
        fmt::print(stderr, "[BioC++ I/O error] {}\n", ext.what());
        return -1;
    }
    catch (decovar_error const & ext)
    {
        fmt::print(stderr, "[deCoVar error] {}\n", ext.what());
        return -1;
    }
#endif
    return 0;
}
