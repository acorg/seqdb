#pragma once

#include "acmacs-base/json-writer.hh"
#include "json-keys.hh"

// ----------------------------------------------------------------------

template <typename Value> class _if_not_empty
{
 public:
    inline _if_not_empty(SeqdbJsonKey key, Value value) : mKey(key), mValue(value) {}

    template <typename RW> friend inline JsonWriterT<RW>& operator <<(JsonWriterT<RW>& writer, const _if_not_empty<Value>& data)
        {
            if (!data.mValue.empty())
                writer << data.mKey << data.mValue;
            return writer;
        }

 private:
    SeqdbJsonKey mKey;
    Value mValue;
};

template <typename Value> _if_not_empty<Value> if_not_empty(SeqdbJsonKey key, Value value) { return _if_not_empty<Value>(key, value); }

// ----------------------------------------------------------------------

template <typename RW> inline JsonWriterT<RW>& operator <<(JsonWriterT<RW>& writer, SeqdbJsonKey key) { const char k = static_cast<char>(key); writer.Key(&k, 1, false); return writer; }

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
