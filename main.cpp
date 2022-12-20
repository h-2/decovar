#include <bio/io/var/reader.hpp>
#include <bio/io/var/writer.hpp>

int main(int argc, char **argv)
{
    bio::io::var::reader_options options{.record = bio::io::var_io::record_default{}};
    bio::io::var::reader reader{argv[1], options};
    bio::io::var::writer writer{argv[2]};


    for (auto & record : reader)
    {
        for (auto && [ id, value ] : record.info)
        {
            if (id == "AF")
            {
                float f = std::get<float>(value);
                if (f > 1e-3) //TODO OPTION
                    writer.push_back(record);
            }
        }
    }
}
