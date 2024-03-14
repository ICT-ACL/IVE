#pragma once

enum OptionKeyword {
    QueryGraphFile = 0,     // -q, The query graph file path
    DataGraphFile = 1,      // -d, The data graph file path
    MaxOutputEmbeddingNum = 10, // -num, The maximum output embedding num
};

class MatchingCommand {
private:
    std::vector<std::string> tokens_;
    std::map<OptionKeyword, std::string> options_key;
    std::map<OptionKeyword, std::string> options_value;

private:
    void processOptions();

public:
    MatchingCommand(int argc, char **argv);

    const std::string getCommandOption(const std::string &option) const {
        std::vector<std::string>::const_iterator itr;
        itr = find(tokens_.begin(), tokens_.end(), option);
        if (itr != tokens_.end() && ++itr != tokens_.end()) {
            return *itr;
        }
        return "";
    }

    bool commandOptionExists(const std::string &option) const {
        return find(tokens_.begin(), tokens_.end(), option) != tokens_.end();
    }

    std::string getDataGraphFilePath() {
        return options_value[OptionKeyword::DataGraphFile];
    }

    std::string getQueryGraphFilePath() {
        return options_value[OptionKeyword::QueryGraphFile];
    }

    std::string getMaximumEmbeddingNum() {
        return options_value[OptionKeyword::MaxOutputEmbeddingNum] == "" ? "MAX" : options_value[OptionKeyword::MaxOutputEmbeddingNum];
    }
};

MatchingCommand::MatchingCommand(int argc, char **argv) {
    for (int i = 1; i < argc; ++i)
        tokens_.push_back(std::string(argv[i]));

    // Initialize options value
    options_key[OptionKeyword::QueryGraphFile] = "-q";
    options_key[OptionKeyword::DataGraphFile] = "-d";
    options_key[OptionKeyword::MaxOutputEmbeddingNum] = "-num";
    processOptions();
};

void MatchingCommand::processOptions() {
    // Query graph file path
    options_value[OptionKeyword::QueryGraphFile] = getCommandOption(options_key[OptionKeyword::QueryGraphFile]);;
    // Data graph file path
    options_value[OptionKeyword::DataGraphFile] = getCommandOption(options_key[OptionKeyword::DataGraphFile]);
    // Maximum output embedding num.
    options_value[OptionKeyword::MaxOutputEmbeddingNum] = getCommandOption(options_key[OptionKeyword::MaxOutputEmbeddingNum]);
}