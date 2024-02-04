using CodecZstd
using JSON

function read_json_zstd(file::String)
    f=open(file) # do f
    stream=ZstdDecompressorStream(f)

    dat = JSON.parse(stream);
    close(stream)
    close(f)
    dat
end