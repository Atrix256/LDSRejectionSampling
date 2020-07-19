#include <random>
#include <vector>
#include <string>
#include <array>

#define COUNT_OF(x) ((sizeof(x)/sizeof(0[x])) / ((size_t)(!(sizeof(x) % sizeof(0[x])))))

#define DETERMINISTIC() true
#define RANDOMIZE_LDS_START() true

static const size_t c_numHistogramBuckets = 10;
static const size_t c_numTests = 1000; // how many times to do the tests, to get average and std dev

// what points in the input stream to dump a report
static const size_t c_reportValues[] =
{
    100,
    1000,
    10000
};
static const size_t c_numReportValues = COUNT_OF(c_reportValues);
static const size_t c_maxReportValue = c_reportValues[c_numReportValues - 1];

// ===================================== Constants =====================================

static const float c_goldenRatio = 1.61803398875f;
static const float c_goldenRatioConjugate = 0.61803398875f;
static const float c_root2 = 1.41421356237f;
static const float c_fractRoot2 = 0.41421356237f;

// generalized golden ratio, from:
// http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/
const float c_R2_g = 1.32471795724474602596f;
const float c_R2_a1 = 1.0f / c_R2_g;
const float c_R2_a2 = 1.0f / (c_R2_g * c_R2_g);

// ===================================== Utils =====================================

typedef std::vector<std::vector<std::string>> CSV;
typedef std::array<float, 2> Vec2;

void SetCSV(CSV& csv, size_t x, size_t y, const char* text)
{
    if (csv.size() <= y)
        csv.resize(y + 1);

    if (csv[y].size() <= x)
        csv[y].resize(x + 1);

    csv[y][x] = text;
}

size_t GetCSVCols(const CSV& csv)
{
    if (csv.size() == 0)
        return 0;
    return csv[0].size();
}

std::mt19937 GetRNG()
{
#if DETERMINISTIC()
    static int seed = 1336;
    seed++;
    std::mt19937 rng(seed);
#else
    std::random_device rd;
    std::mt19937 rng(rd());
#endif
    return rng;
}

float fract(float f)
{
    return f - floor(f);
}

int FloatToItem(float f, int numItems)
{
    // if f is in [0,1], remaps to [0, numItems-1]
    return std::min(int(f * float(numItems)), numItems - 1);
}

float Lerp(float a, float b, float t)
{
    return a * (1.0f - t) + b * t;
}

// ===================================== LDSs =====================================

static size_t Ruler(size_t n)
{
    size_t ret = 0;
    while (n != 0 && (n & 1) == 0)
    {
        n /= 2;
        ++ret;
    }
    return ret;
}

void Sobol(std::vector<Vec2>& values, size_t numValues)
{
    // x axis
    values.resize(numValues);
    size_t sampleInt = 0;
    for (size_t i = 0; i < numValues; ++i)
    {
        size_t ruler = Ruler(i + 1);
        size_t direction = size_t(size_t(1) << size_t(31 - ruler));
        sampleInt = sampleInt ^ direction;
        values[i][0] = float(sampleInt) / std::pow(2.0f, 32.0f);
    }

    // y axis
    // Code adapted from http://web.maths.unsw.edu.au/~fkuo/sobol/
    // uses numbers: new-joe-kuo-6.21201

    // Direction numbers
    std::vector<size_t> V;
    V.resize((size_t)ceil(log((double)numValues + 1) / log(2.0)));  //+1 because we are skipping index 0
    V[0] = size_t(1) << size_t(31);
    for (size_t i = 1; i < V.size(); ++i)
        V[i] = V[i - 1] ^ (V[i - 1] >> 1);

    // Samples
    sampleInt = 0;
    for (size_t i = 0; i < numValues; ++i) {
        size_t ruler = Ruler(i + 1);
        sampleInt = sampleInt ^ V[ruler];
        values[i][1] = float(sampleInt) / std::pow(2.0f, 32.0f);
    }
}

// ===================================== PDFs =====================================

namespace PDF
{
    // PDF: y = 1
    struct UniformWhite
    {
        UniformWhite()
        {
            m_rng = GetRNG();
        }

        static float PF(float x)
        {
            return 1.0f;
        }

        float Generate(float generatorPFMultiplier)
        {
            m_attempts++;
            std::uniform_real_distribution<float> dist(0.0f, 1.0f);
            return dist(m_rng);
        }

        size_t m_attempts = 0;
        std::mt19937 m_rng;
    };

    // PDF: y = 1
    struct UniformLDS_GR
    {
        UniformLDS_GR()
        {
            #if RANDOMIZE_LDS_START()
            std::mt19937 rng = GetRNG();
            std::uniform_real_distribution<float> dist(0.0f, 1.0f);
            m_value = dist(rng);
            #else
            m_value = 0.0f;
            #endif
        }

        static float PF(float x)
        {
            return 1.0f;
        }

        float Generate(float generatorPFMultiplier)
        {
            m_attempts++;
            m_value = fract(m_value + c_goldenRatioConjugate);
            return m_value;
        }

        size_t m_attempts = 0;
        float m_value;
    };

    // PDF: y = 1
    struct UniformLDS_Root2
    {
        UniformLDS_Root2()
        {
            #if RANDOMIZE_LDS_START()
            std::mt19937 rng = GetRNG();
            std::uniform_real_distribution<float> dist(0.0f, 1.0f);
            m_value = dist(rng);
            #else
            m_value = 0.0f;
            #endif
        }

        static float PF(float x)
        {
            return 1.0f;
        }

        float Generate(float generatorPFMultiplier)
        {
            m_attempts++;
            m_value = fract(m_value + c_fractRoot2);
            return m_value;
        }

        size_t m_attempts = 0;
        float m_value;
    };

    // PDF: y = 1
    struct UniformLDS_R2_x
    {
        UniformLDS_R2_x()
        {
            #if RANDOMIZE_LDS_START()
            std::mt19937 rng = GetRNG();
            std::uniform_real_distribution<float> dist(0.0f, 1.0f);
            m_value = dist(rng);
            #else
            m_value = 0.0f;
            #endif
        }

        static float PF(float x)
        {
            return 1.0f;
        }

        float Generate(float generatorPFMultiplier)
        {
            m_attempts++;
            m_value = fract(m_value + c_R2_a1);
            return m_value;
        }

        size_t m_attempts = 0;
        float m_value;
    };

    // PDF: y = 1
    struct UniformLDS_R2_y
    {
        UniformLDS_R2_y()
        {
            #if RANDOMIZE_LDS_START()
            std::mt19937 rng = GetRNG();
            std::uniform_real_distribution<float> dist(0.0f, 1.0f);
            m_value = dist(rng);
            #else
            m_value = 0.0f;
            #endif
        }

        static float PF(float x)
        {
            return 1.0f;
        }

        float Generate(float generatorPFMultiplier)
        {
            m_attempts++;
            m_value = fract(m_value + c_R2_a2);
            return m_value;
        }

        size_t m_attempts = 0;
        float m_value;
    };

    // PDF: y = 1
    struct UniformLDS_Sobol_x
    {
        UniformLDS_Sobol_x()
        {
            #if RANDOMIZE_LDS_START()
            std::mt19937 rng = GetRNG();
            std::uniform_real_distribution<float> dist(0.0f, 1.0f);
            m_rngOffset = dist(rng);
            #else
            m_rngOffset = 0.0f;
            #endif
        }

        static float PF(float x)
        {
            return 1.0f;
        }

        float Generate(float generatorPFMultiplier)
        {
            if (m_index >= m_sobolPoints.size())
                Sobol(m_sobolPoints, m_sobolPoints.size() + 1000);

            m_attempts++;
            float ret = m_sobolPoints[m_index][0];
            m_index++;
            return fract(ret + m_rngOffset);
        }

        size_t m_index = 0;
        size_t m_attempts = 0;
        std::vector<Vec2> m_sobolPoints;
        float m_rngOffset;
    };

    // PDF: y = 1
    struct UniformLDS_Sobol_y
    {
        UniformLDS_Sobol_y()
        {
            #if RANDOMIZE_LDS_START()
            std::mt19937 rng = GetRNG();
            std::uniform_real_distribution<float> dist(0.0f, 1.0f);
            m_rngOffset = dist(rng);
            #else
            m_rngOffset = 0.0f;
            #endif
        }

        static float PF(float x)
        {
            return 1.0f;
        }

        float Generate(float generatorPFMultiplier)
        {
            if (m_index >= m_sobolPoints.size())
                Sobol(m_sobolPoints, m_sobolPoints.size() + 1000);

            m_attempts++;
            float ret = m_sobolPoints[m_index][1];
            m_index++;
            return fract(ret + m_rngOffset);
        }

        size_t m_index = 0;
        size_t m_attempts = 0;
        std::vector<Vec2> m_sobolPoints;
        float m_rngOffset;
    };


    // PDF:  y = (2x + 3) / 4
    // CDF:  y = (x^2 + 3x) / 4 
    // PF:   y = (2x + 3) / 5
    template <typename Generator, typename Validator>
    struct Linear
    {
        static float CDF(float x)
        {
            return (x*x + 3.0f*x) / 4.0f;
        }

        static float PF(float x)
        {
            return (2.0f * x + 3.0f) / 5.0f;
        }

        float Generate(float generatorPFMultiplier)
        {
            while (true)
            {
                m_attempts++;
                float x = m_generator.Generate(1.0f); // doesn't handle nested non-uniform PDFs but whatever
                float probability = PF(x) / (Generator::PF(x) * generatorPFMultiplier);
                if (m_validator.Generate(1.0f) <= probability)
                    return x;
            }
        }

        size_t m_attempts = 0;
        Generator m_generator;
        Validator m_validator;
    };

    // PDF:  y = (x^3 -10x^2 + 5x + 11) / 10.417
    // CDF:  y = 0.0959969  * (11x + 2.5x^2 - 3.33333x^3 + 0.25x^4)
    // PF:   y = (x^3 -10x^2 + 5x + 11) / 12
    template <typename Generator, typename Validator>
    struct Cubic
    {
        static float CDF(float x)
        {
            return 0.0959969f * (11.0f * x + 2.5f*x*x - 3.33333f*x*x*x + 0.25f * x*x*x*x);
        }

        static float PF(float x)
        {
            return (x*x*x - 10.0f * x*x + 5.0f * x + 11.0f) / 12.0f;
        }

        float Generate(float generatorPFMultiplier)
        {
            while (true)
            {
                m_attempts++;
                float x = m_generator.Generate(1.0f); // doesn't handle nested non-uniform PDFs but whatever
                float probability = PF(x) / (Generator::PF(x) * generatorPFMultiplier);
                if (m_validator.Generate(1.0f) <= probability)
                    return x;
            }
        }

        size_t m_attempts = 0;
        Generator m_generator;
        Validator m_validator;
    };
};

// ===================================== Code =====================================

template <typename TPDF>
void GenerateSequence(std::vector<float>& sequence, size_t count, float generatorPFMultiplier, std::vector<size_t>& attempts)
{
    TPDF pdf;
    sequence.resize(count);
    attempts.resize(0);
    for (float& f : sequence)
    {
        f = pdf.Generate(generatorPFMultiplier);
        attempts.push_back(pdf.m_attempts);
    }
}

void CalculateHistogram(const std::vector<float>& samples, size_t sampleCount, std::vector<float>& histogram)
{
    histogram.resize(c_numHistogramBuckets);
    std::fill(histogram.begin(), histogram.end(), 0.0f);
    for (size_t index = 0; index < sampleCount; ++index)
        histogram[FloatToItem(samples[index], c_numHistogramBuckets)] += 1.0f / float(sampleCount);
}

void WriteHistogram(FILE* file, const char* label, const std::vector<float>& histogram)
{
    fprintf(file, "\"%s\"", label);
    for (float f : histogram)
        fprintf(file, ",\"%f\"", f);
    fprintf(file, "\n");
}

void WriteHistogramError(FILE* file, const char* label, const std::vector<float>& histogram, const std::vector<float>& expectedHistogram)
{
    fprintf(file, "\"%s\"", label);
    for (size_t i = 0; i < c_numHistogramBuckets; ++i)
        fprintf(file, ",\"%f\"", histogram[i] - expectedHistogram[i]);
    fprintf(file, "\n");
}

void WriteHistogramStdDev(FILE* file, const char* label, const std::vector<float>& histogramAvg, const std::vector<float>& histogramSqAvg)
{
    fprintf(file, "\"%s\"", label);
    for (size_t i = 0; i < c_numHistogramBuckets; ++i)
        fprintf(file, ",\"%f\"", sqrt(abs(histogramSqAvg[i] - histogramAvg[i] * histogramAvg[i])));
    fprintf(file, "\n");
}

void WriteSurvival(CSV& csv, const char* label, const std::vector<size_t>& attempts, size_t count)
{
    size_t col = GetCSVCols(csv);
    SetCSV(csv, col, 0, label);
    char buffer[256];
    for (size_t index = 0; index < count; ++index)
    {
        sprintf_s(buffer, "%zu", attempts[index]);
        SetCSV(csv, col, index + 1, buffer);
    }
}

void WriteSurvival(CSV& csv, const char* label, const std::vector<float>& attempts, size_t count)
{
    size_t col = GetCSVCols(csv);
    SetCSV(csv, col, 0, label);
    char buffer[256];
    for (size_t index = 0; index < count; ++index)
    {
        sprintf_s(buffer, "%f", attempts[index]);
        SetCSV(csv, col, index + 1, buffer);
    }
}

void WriteSurvivalStdDev(CSV& csv, const char* label, const std::vector<float>& attemptsAvg, const std::vector<float>& attemptsSqAvg, size_t count)
{
    size_t col = GetCSVCols(csv);
    SetCSV(csv, col, 0, label);
    char buffer[256];
    for (size_t index = 0; index < count; ++index)
    {
        float avg = attemptsAvg[index];
        float sqAvg = attemptsSqAvg[index];
        float stdDev = sqrt(abs(sqAvg - avg * avg));
        sprintf_s(buffer, "%f", stdDev);
        SetCSV(csv, col, index + 1, buffer);
    }
}

void WriteCSV(const CSV& csv, const char* fileName)
{
    FILE* file = nullptr;
    fopen_s(&file, fileName, "w+t");
    for (const std::vector<std::string>& row : csv)
    {
        bool first = true;
        for (const std::string& item : row)
        {
            fprintf(file, "%s\"%s\"", first ? "" : ",", item.c_str());
            first = false;
        }
        fprintf(file, "\n");
    }
    fclose(file);
}

void CombineAttemptsData(const std::vector<size_t>& attempts, std::vector<float>& attemptsAvg, std::vector<float>& attemptsSqAvg, size_t sampleIndex)
{
    attemptsAvg.resize(attempts.size());
    attemptsSqAvg.resize(attempts.size());
    for (size_t index = 0; index < attempts.size(); ++index)
    {
        attemptsAvg[index] = Lerp(attemptsAvg[index], float(attempts[index]), 1.0f / float(sampleIndex + 1));
        attemptsSqAvg[index] = Lerp(attemptsSqAvg[index], float(attempts[index])*float(attempts[index]), 1.0f / float(sampleIndex + 1));
    }
}

void HistogramCombine(const std::vector<float>& histogram, std::vector<float>& histogramAvg, std::vector<float>& histogramSqAvg, size_t sampleIndex)
{
    if (sampleIndex == 0)
    {
        histogramAvg.resize(c_numHistogramBuckets, 0.0f);
        histogramSqAvg.resize(c_numHistogramBuckets, 0.0f);
    }

    for (size_t index = 0; index < c_numHistogramBuckets; ++index)
    {
        float value = histogram[index];
        histogramAvg[index] = Lerp(histogramAvg[index], value, 1.0f / float(sampleIndex + 1));
        histogramSqAvg[index] = Lerp(histogramSqAvg[index], value*value, 1.0f / float(sampleIndex + 1));
    }
}

int main(int argc, char ** argv)
{
    // uniform to linear
    {
        // calculate the expected histogram
        std::vector<float> expectedHistogram;
        for (size_t bucketIndex = 0; bucketIndex < c_numHistogramBuckets; ++bucketIndex)
        {
            float p0 = float(bucketIndex) / float(c_numHistogramBuckets);
            float p1 = float(bucketIndex + 1) / float(c_numHistogramBuckets);

            float cdf0 = PDF::Linear<PDF::UniformWhite, PDF::UniformWhite>::CDF(p0);
            float cdf1 = PDF::Linear<PDF::UniformWhite, PDF::UniformWhite>::CDF(p1);

            expectedHistogram.push_back(cdf1 - cdf0);
        }

        std::vector<float> histogram[c_numReportValues][6];
        std::vector<float> histogramAvg[c_numReportValues][6];
        std::vector<float> histogramSqAvg[c_numReportValues][6];

        std::vector<size_t> attempts[6];
        std::vector<float> attemptsAvg[6];
        std::vector<float> attemptsSqAvg[6];

        for (size_t testIndex = 0; testIndex < c_numTests; ++testIndex)
        {
            // generate the samples
            std::vector<float> samples[6];
            GenerateSequence<PDF::Linear<PDF::UniformWhite, PDF::UniformWhite>>(samples[0], c_maxReportValue, 1.0f, attempts[0]);
            GenerateSequence<PDF::Linear<PDF::UniformWhite, PDF::UniformLDS_GR>>(samples[1], c_maxReportValue, 1.0f, attempts[1]);
            GenerateSequence<PDF::Linear<PDF::UniformLDS_Root2, PDF::UniformWhite>>(samples[2], c_maxReportValue, 1.0f, attempts[2]);
            GenerateSequence<PDF::Linear<PDF::UniformLDS_Root2, PDF::UniformLDS_GR>>(samples[3], c_maxReportValue, 1.0f, attempts[3]);
            GenerateSequence<PDF::Linear<PDF::UniformLDS_R2_x, PDF::UniformLDS_R2_y>>(samples[4], c_maxReportValue, 1.0f, attempts[4]);
            GenerateSequence<PDF::Linear<PDF::UniformLDS_Sobol_x, PDF::UniformLDS_Sobol_y>>(samples[5], c_maxReportValue, 1.0f, attempts[5]);

            // combine the attempts data
            CombineAttemptsData(attempts[0], attemptsAvg[0], attemptsSqAvg[0], testIndex);
            CombineAttemptsData(attempts[1], attemptsAvg[1], attemptsSqAvg[1], testIndex);
            CombineAttemptsData(attempts[2], attemptsAvg[2], attemptsSqAvg[2], testIndex);
            CombineAttemptsData(attempts[3], attemptsAvg[3], attemptsSqAvg[3], testIndex);
            CombineAttemptsData(attempts[4], attemptsAvg[4], attemptsSqAvg[4], testIndex);
            CombineAttemptsData(attempts[5], attemptsAvg[5], attemptsSqAvg[5], testIndex);

            // calculate data for each report
            for (size_t reportIndex = 0; reportIndex < c_numReportValues; ++reportIndex)
            {
                size_t sampleCount = c_reportValues[reportIndex];

                // calculate the histograms
                CalculateHistogram(samples[0], sampleCount, histogram[reportIndex][0]);
                CalculateHistogram(samples[1], sampleCount, histogram[reportIndex][1]);
                CalculateHistogram(samples[2], sampleCount, histogram[reportIndex][2]);
                CalculateHistogram(samples[3], sampleCount, histogram[reportIndex][3]);
                CalculateHistogram(samples[4], sampleCount, histogram[reportIndex][4]);
                CalculateHistogram(samples[5], sampleCount, histogram[reportIndex][5]);

                // combine the histograms
                HistogramCombine(histogram[reportIndex][0], histogramAvg[reportIndex][0], histogramSqAvg[reportIndex][0], testIndex);
                HistogramCombine(histogram[reportIndex][1], histogramAvg[reportIndex][1], histogramSqAvg[reportIndex][1], testIndex);
                HistogramCombine(histogram[reportIndex][2], histogramAvg[reportIndex][2], histogramSqAvg[reportIndex][2], testIndex);
                HistogramCombine(histogram[reportIndex][3], histogramAvg[reportIndex][3], histogramSqAvg[reportIndex][3], testIndex);
                HistogramCombine(histogram[reportIndex][4], histogramAvg[reportIndex][4], histogramSqAvg[reportIndex][4], testIndex);
                HistogramCombine(histogram[reportIndex][5], histogramAvg[reportIndex][5], histogramSqAvg[reportIndex][5], testIndex);
            }
        }

        // write the attempts, to show the efficiency of the rejection sampling
        {
            CSV csv;
            SetCSV(csv, GetCSVCols(csv), 0, "Survival");
            WriteSurvival(csv, "white/white", attempts[0], 100);
            WriteSurvival(csv, "white/LDS", attempts[1], 100);
            WriteSurvival(csv, "LDS/white", attempts[2], 100);
            WriteSurvival(csv, "LDS/LDS", attempts[3], 100);
            WriteSurvival(csv, "R2 LDS", attempts[4], 100);
            WriteSurvival(csv, "Sobol LDS", attempts[5], 100);

            SetCSV(csv, GetCSVCols(csv), 0, "Survival Avg");
            WriteSurvival(csv, "white/white", attemptsAvg[0], 100);
            WriteSurvival(csv, "white/LDS", attemptsAvg[1], 100);
            WriteSurvival(csv, "LDS/white", attemptsAvg[2], 100);
            WriteSurvival(csv, "LDS/LDS", attemptsAvg[3], 100);
            WriteSurvival(csv, "R2 LDS", attemptsAvg[4], 100);
            WriteSurvival(csv, "Sobol LDS", attemptsAvg[5], 100);

            SetCSV(csv, GetCSVCols(csv), 0, "Survival Std Dev");
            WriteSurvivalStdDev(csv, "white/white", attemptsAvg[0], attemptsSqAvg[0], 100);
            WriteSurvivalStdDev(csv, "white/LDS", attemptsAvg[1], attemptsSqAvg[1], 100);
            WriteSurvivalStdDev(csv, "LDS/white", attemptsAvg[2], attemptsSqAvg[2], 100);
            WriteSurvivalStdDev(csv, "LDS/LDS", attemptsAvg[3], attemptsSqAvg[3], 100);
            WriteSurvivalStdDev(csv, "R2 LDS", attemptsAvg[4], attemptsSqAvg[4], 100);
            WriteSurvivalStdDev(csv, "Sobol LDS", attemptsAvg[5], attemptsSqAvg[5], 100);

            WriteCSV(csv, "out/uni_lin_survival.csv");
        }

        // make the reports
        for (size_t reportIndex = 0; reportIndex < c_numReportValues; ++reportIndex)
        {
            size_t sampleCount = c_reportValues[reportIndex];

            // open the file
            char buffer[256];
            sprintf_s(buffer, "out/uni_lin_%zu.csv", sampleCount);
            FILE* file = nullptr;
            fopen_s(&file, buffer, "w+t");

            // write the expected histogram
            fprintf(file, "\"Expected\"");
            for (float f : expectedHistogram)
                fprintf(file, ",\"%f\"", f);
            fprintf(file, "\n");

            // write histograms
            WriteHistogram(file, "white/white", histogram[reportIndex][0]);
            WriteHistogram(file, "white/LDS", histogram[reportIndex][1]);
            WriteHistogram(file, "LDS/white", histogram[reportIndex][2]);
            WriteHistogram(file, "LDS/LDS", histogram[reportIndex][3]);
            WriteHistogram(file, "R2 LDS", histogram[reportIndex][4]);
            WriteHistogram(file, "Sobol LDS", histogram[reportIndex][5]);

            // write error
            fprintf(file, "\n\"Error:\"\n");
            WriteHistogramError(file, "white/white", histogram[reportIndex][0], expectedHistogram);
            WriteHistogramError(file, "white/LDS", histogram[reportIndex][1], expectedHistogram);
            WriteHistogramError(file, "LDS/white", histogram[reportIndex][2], expectedHistogram);
            WriteHistogramError(file, "LDS/LDS", histogram[reportIndex][3], expectedHistogram);
            WriteHistogramError(file, "R2 LDS", histogram[reportIndex][4], expectedHistogram);
            WriteHistogramError(file, "Sobol LDS", histogram[reportIndex][5], expectedHistogram);

            // write average error
            fprintf(file, "\n\"Avg Error:\"\n");
            WriteHistogramError(file, "white/white", histogramAvg[reportIndex][0], expectedHistogram);
            WriteHistogramError(file, "white/LDS", histogramAvg[reportIndex][1], expectedHistogram);
            WriteHistogramError(file, "LDS/white", histogramAvg[reportIndex][2], expectedHistogram);
            WriteHistogramError(file, "LDS/LDS", histogramAvg[reportIndex][3], expectedHistogram);
            WriteHistogramError(file, "R2 LDS", histogramAvg[reportIndex][4], expectedHistogram);
            WriteHistogramError(file, "Sobol LDS", histogramAvg[reportIndex][5], expectedHistogram);

            // write std dev
            fprintf(file, "\n\"Std Dev:\"\n");
            WriteHistogramStdDev(file, "white/white", histogramAvg[reportIndex][0], histogramSqAvg[reportIndex][0]);
            WriteHistogramStdDev(file, "white/LDS", histogramAvg[reportIndex][1], histogramSqAvg[reportIndex][1]);
            WriteHistogramStdDev(file, "LDS/white", histogramAvg[reportIndex][2], histogramSqAvg[reportIndex][2]);
            WriteHistogramStdDev(file, "LDS/LDS", histogramAvg[reportIndex][3], histogramSqAvg[reportIndex][3]);
            WriteHistogramStdDev(file, "R2 LDS", histogramAvg[reportIndex][4], histogramSqAvg[reportIndex][4]);
            WriteHistogramStdDev(file, "Sobol LDS", histogramAvg[reportIndex][5], histogramSqAvg[reportIndex][5]);

            fclose(file);
        }
    }

    // [uniform to linear] to cubic
    {
        // calculate the expected histogram
        std::vector<float> expectedHistogram;
        for (size_t bucketIndex = 0; bucketIndex < c_numHistogramBuckets; ++bucketIndex)
        {
            float p0 = float(bucketIndex) / float(c_numHistogramBuckets);
            float p1 = float(bucketIndex + 1) / float(c_numHistogramBuckets);

            float cdf0 = PDF::Cubic<PDF::Linear<PDF::UniformWhite, PDF::UniformWhite>, PDF::UniformWhite>::CDF(p0);
            float cdf1 = PDF::Cubic<PDF::Linear<PDF::UniformWhite, PDF::UniformWhite>, PDF::UniformWhite>::CDF(p1);

            expectedHistogram.push_back(cdf1 - cdf0);
        }

        std::vector<float> histogram[c_numReportValues][6];
        std::vector<float> histogramAvg[c_numReportValues][6];
        std::vector<float> histogramSqAvg[c_numReportValues][6];

        std::vector<size_t> attempts[6];
        std::vector<float> attemptsAvg[6];
        std::vector<float> attemptsSqAvg[6];

        for (size_t testIndex = 0; testIndex < c_numTests; ++testIndex)
        {
            // generate the samples
            std::vector<float> samples[6];
            GenerateSequence<PDF::Cubic<PDF::Linear<PDF::UniformWhite, PDF::UniformWhite>, PDF::UniformWhite>>(samples[0], c_maxReportValue, 1.6f, attempts[0]);
            GenerateSequence<PDF::Cubic<PDF::Linear<PDF::UniformWhite, PDF::UniformLDS_GR>, PDF::UniformLDS_GR>>(samples[1], c_maxReportValue, 1.6f, attempts[1]);
            GenerateSequence<PDF::Cubic<PDF::Linear<PDF::UniformLDS_Root2, PDF::UniformWhite>, PDF::UniformWhite>>(samples[2], c_maxReportValue, 1.6f, attempts[2]);
            GenerateSequence<PDF::Cubic<PDF::Linear<PDF::UniformLDS_Root2, PDF::UniformLDS_GR>, PDF::UniformLDS_GR>>(samples[3], c_maxReportValue, 1.6f, attempts[3]);
            GenerateSequence<PDF::Cubic<PDF::Linear<PDF::UniformLDS_R2_x, PDF::UniformLDS_R2_y>, PDF::UniformLDS_R2_y>>(samples[4], c_maxReportValue, 1.6f, attempts[4]);
            GenerateSequence<PDF::Cubic<PDF::Linear<PDF::UniformLDS_Sobol_x, PDF::UniformLDS_Sobol_y>, PDF::UniformLDS_Sobol_y>>(samples[5], c_maxReportValue, 1.6f, attempts[5]);

            // combine the attempts data
            CombineAttemptsData(attempts[0], attemptsAvg[0], attemptsSqAvg[0], testIndex);
            CombineAttemptsData(attempts[1], attemptsAvg[1], attemptsSqAvg[1], testIndex);
            CombineAttemptsData(attempts[2], attemptsAvg[2], attemptsSqAvg[2], testIndex);
            CombineAttemptsData(attempts[3], attemptsAvg[3], attemptsSqAvg[3], testIndex);
            CombineAttemptsData(attempts[4], attemptsAvg[4], attemptsSqAvg[4], testIndex);
            CombineAttemptsData(attempts[5], attemptsAvg[5], attemptsSqAvg[5], testIndex);

            // calculate data for each report
            for (size_t reportIndex = 0; reportIndex < c_numReportValues; ++reportIndex)
            {
                size_t sampleCount = c_reportValues[reportIndex];

                // calculate the histograms
                CalculateHistogram(samples[0], sampleCount, histogram[reportIndex][0]);
                CalculateHistogram(samples[1], sampleCount, histogram[reportIndex][1]);
                CalculateHistogram(samples[2], sampleCount, histogram[reportIndex][2]);
                CalculateHistogram(samples[3], sampleCount, histogram[reportIndex][3]);
                CalculateHistogram(samples[4], sampleCount, histogram[reportIndex][4]);
                CalculateHistogram(samples[5], sampleCount, histogram[reportIndex][5]);

                // combine the histograms
                HistogramCombine(histogram[reportIndex][0], histogramAvg[reportIndex][0], histogramSqAvg[reportIndex][0], testIndex);
                HistogramCombine(histogram[reportIndex][1], histogramAvg[reportIndex][1], histogramSqAvg[reportIndex][1], testIndex);
                HistogramCombine(histogram[reportIndex][2], histogramAvg[reportIndex][2], histogramSqAvg[reportIndex][2], testIndex);
                HistogramCombine(histogram[reportIndex][3], histogramAvg[reportIndex][3], histogramSqAvg[reportIndex][3], testIndex);
                HistogramCombine(histogram[reportIndex][4], histogramAvg[reportIndex][4], histogramSqAvg[reportIndex][4], testIndex);
                HistogramCombine(histogram[reportIndex][5], histogramAvg[reportIndex][5], histogramSqAvg[reportIndex][5], testIndex);
            }
        }

        // write the attempts, to show the efficiency of the rejection sampling
        {
            CSV csv;
            SetCSV(csv, GetCSVCols(csv), 0, "Survival");
            WriteSurvival(csv, "white/white", attempts[0], 100);
            WriteSurvival(csv, "white/LDS", attempts[1], 100);
            WriteSurvival(csv, "LDS/white", attempts[2], 100);
            WriteSurvival(csv, "LDS/LDS", attempts[3], 100);
            WriteSurvival(csv, "R2 LDS", attempts[4], 100);
            WriteSurvival(csv, "Sobol LDS", attempts[5], 100);

            SetCSV(csv, GetCSVCols(csv), 0, "Survival Avg");
            WriteSurvival(csv, "white/white", attemptsAvg[0], 100);
            WriteSurvival(csv, "white/LDS", attemptsAvg[1], 100);
            WriteSurvival(csv, "LDS/white", attemptsAvg[2], 100);
            WriteSurvival(csv, "LDS/LDS", attemptsAvg[3], 100);
            WriteSurvival(csv, "R2 LDS", attemptsAvg[4], 100);
            WriteSurvival(csv, "Sobol LDS", attemptsAvg[5], 100);

            SetCSV(csv, GetCSVCols(csv), 0, "Survival Std Dev");
            WriteSurvivalStdDev(csv, "white/white", attemptsAvg[0], attemptsSqAvg[0], 100);
            WriteSurvivalStdDev(csv, "white/LDS", attemptsAvg[1], attemptsSqAvg[1], 100);
            WriteSurvivalStdDev(csv, "LDS/white", attemptsAvg[2], attemptsSqAvg[2], 100);
            WriteSurvivalStdDev(csv, "LDS/LDS", attemptsAvg[3], attemptsSqAvg[3], 100);
            WriteSurvivalStdDev(csv, "R2 LDS", attemptsAvg[4], attemptsSqAvg[4], 100);
            WriteSurvivalStdDev(csv, "Sobol LDS", attemptsAvg[5], attemptsSqAvg[5], 100);

            WriteCSV(csv, "out/lin_cub_survival.csv");
        }

        // make the reports
        for (size_t reportIndex = 0; reportIndex < c_numReportValues; ++reportIndex)
        {
            size_t sampleCount = c_reportValues[reportIndex];

            // open the file
            char buffer[256];
            sprintf_s(buffer, "out/lin_cub_%zu.csv", sampleCount);
            FILE* file = nullptr;
            fopen_s(&file, buffer, "w+t");

            // write the expected histogram
            fprintf(file, "\"Expected\"");
            for (float f : expectedHistogram)
                fprintf(file, ",\"%f\"", f);
            fprintf(file, "\n");

            // write histograms
            WriteHistogram(file, "white/white", histogram[reportIndex][0]);
            WriteHistogram(file, "white/LDS", histogram[reportIndex][1]);
            WriteHistogram(file, "LDS/white", histogram[reportIndex][2]);
            WriteHistogram(file, "LDS/LDS", histogram[reportIndex][3]);
            WriteHistogram(file, "R2 LDS", histogram[reportIndex][4]);
            WriteHistogram(file, "Sobol LDS", histogram[reportIndex][5]);

            // write error
            fprintf(file, "\n\"Error:\"\n");
            WriteHistogramError(file, "white/white", histogram[reportIndex][0], expectedHistogram);
            WriteHistogramError(file, "white/LDS", histogram[reportIndex][1], expectedHistogram);
            WriteHistogramError(file, "LDS/white", histogram[reportIndex][2], expectedHistogram);
            WriteHistogramError(file, "LDS/LDS", histogram[reportIndex][3], expectedHistogram);
            WriteHistogramError(file, "R2 LDS", histogram[reportIndex][4], expectedHistogram);
            WriteHistogramError(file, "Sobol LDS", histogram[reportIndex][5], expectedHistogram);

            // write average error
            fprintf(file, "\n\"Avg Error:\"\n");
            WriteHistogramError(file, "white/white", histogramAvg[reportIndex][0], expectedHistogram);
            WriteHistogramError(file, "white/LDS", histogramAvg[reportIndex][1], expectedHistogram);
            WriteHistogramError(file, "LDS/white", histogramAvg[reportIndex][2], expectedHistogram);
            WriteHistogramError(file, "LDS/LDS", histogramAvg[reportIndex][3], expectedHistogram);
            WriteHistogramError(file, "R2 LDS", histogramAvg[reportIndex][4], expectedHistogram);
            WriteHistogramError(file, "Sobol LDS", histogramAvg[reportIndex][5], expectedHistogram);

            // write std dev
            fprintf(file, "\n\"Std Dev:\"\n");
            WriteHistogramStdDev(file, "white/white", histogramAvg[reportIndex][0], histogramSqAvg[reportIndex][0]);
            WriteHistogramStdDev(file, "white/LDS", histogramAvg[reportIndex][1], histogramSqAvg[reportIndex][1]);
            WriteHistogramStdDev(file, "LDS/white", histogramAvg[reportIndex][2], histogramSqAvg[reportIndex][2]);
            WriteHistogramStdDev(file, "LDS/LDS", histogramAvg[reportIndex][3], histogramSqAvg[reportIndex][3]);
            WriteHistogramStdDev(file, "R2 LDS", histogramAvg[reportIndex][4], histogramSqAvg[reportIndex][4]);
            WriteHistogramStdDev(file, "Sobol LDS", histogramAvg[reportIndex][5], histogramSqAvg[reportIndex][5]);

            fclose(file);
        }
    }

    // uniform to cubic
    {
        // calculate the expected histogram
        std::vector<float> expectedHistogram;
        for (size_t bucketIndex = 0; bucketIndex < c_numHistogramBuckets; ++bucketIndex)
        {
            float p0 = float(bucketIndex) / float(c_numHistogramBuckets);
            float p1 = float(bucketIndex + 1) / float(c_numHistogramBuckets);

            float cdf0 = PDF::Cubic<PDF::UniformWhite, PDF::UniformWhite>::CDF(p0);
            float cdf1 = PDF::Cubic<PDF::UniformWhite, PDF::UniformWhite>::CDF(p1);

            expectedHistogram.push_back(cdf1 - cdf0);
        }

        std::vector<float> histogram[c_numReportValues][6];
        std::vector<float> histogramAvg[c_numReportValues][6];
        std::vector<float> histogramSqAvg[c_numReportValues][6];

        std::vector<size_t> attempts[6];
        std::vector<float> attemptsAvg[6];
        std::vector<float> attemptsSqAvg[6];

        for (size_t testIndex = 0; testIndex < c_numTests; ++testIndex)
        {
            // generate the samples
            std::vector<float> samples[6];
            GenerateSequence<PDF::Cubic<PDF::UniformWhite, PDF::UniformWhite>>(samples[0], c_maxReportValue, 1.0f, attempts[0]);
            GenerateSequence<PDF::Cubic<PDF::UniformWhite, PDF::UniformLDS_GR>>(samples[1], c_maxReportValue, 1.0f, attempts[1]);
            GenerateSequence<PDF::Cubic<PDF::UniformLDS_Root2, PDF::UniformWhite>>(samples[2], c_maxReportValue, 1.0f, attempts[2]);
            GenerateSequence<PDF::Cubic<PDF::UniformLDS_Root2, PDF::UniformLDS_GR>>(samples[3], c_maxReportValue, 1.0f, attempts[3]);
            GenerateSequence<PDF::Cubic<PDF::UniformLDS_R2_x, PDF::UniformLDS_R2_y>>(samples[4], c_maxReportValue, 1.0f, attempts[4]);
            GenerateSequence<PDF::Cubic<PDF::UniformLDS_Sobol_x, PDF::UniformLDS_Sobol_y>>(samples[5], c_maxReportValue, 1.0f, attempts[5]);

            // combine the attempts data
            CombineAttemptsData(attempts[0], attemptsAvg[0], attemptsSqAvg[0], testIndex);
            CombineAttemptsData(attempts[1], attemptsAvg[1], attemptsSqAvg[1], testIndex);
            CombineAttemptsData(attempts[2], attemptsAvg[2], attemptsSqAvg[2], testIndex);
            CombineAttemptsData(attempts[3], attemptsAvg[3], attemptsSqAvg[3], testIndex);
            CombineAttemptsData(attempts[4], attemptsAvg[4], attemptsSqAvg[4], testIndex);
            CombineAttemptsData(attempts[5], attemptsAvg[5], attemptsSqAvg[5], testIndex);

            // calculate data for each report
            for (size_t reportIndex = 0; reportIndex < c_numReportValues; ++reportIndex)
            {
                size_t sampleCount = c_reportValues[reportIndex];

                // calculate the histograms
                CalculateHistogram(samples[0], sampleCount, histogram[reportIndex][0]);
                CalculateHistogram(samples[1], sampleCount, histogram[reportIndex][1]);
                CalculateHistogram(samples[2], sampleCount, histogram[reportIndex][2]);
                CalculateHistogram(samples[3], sampleCount, histogram[reportIndex][3]);
                CalculateHistogram(samples[4], sampleCount, histogram[reportIndex][4]);
                CalculateHistogram(samples[5], sampleCount, histogram[reportIndex][5]);

                // combine the histograms
                HistogramCombine(histogram[reportIndex][0], histogramAvg[reportIndex][0], histogramSqAvg[reportIndex][0], testIndex);
                HistogramCombine(histogram[reportIndex][1], histogramAvg[reportIndex][1], histogramSqAvg[reportIndex][1], testIndex);
                HistogramCombine(histogram[reportIndex][2], histogramAvg[reportIndex][2], histogramSqAvg[reportIndex][2], testIndex);
                HistogramCombine(histogram[reportIndex][3], histogramAvg[reportIndex][3], histogramSqAvg[reportIndex][3], testIndex);
                HistogramCombine(histogram[reportIndex][4], histogramAvg[reportIndex][4], histogramSqAvg[reportIndex][4], testIndex);
                HistogramCombine(histogram[reportIndex][5], histogramAvg[reportIndex][5], histogramSqAvg[reportIndex][5], testIndex);
            }
        }

        // write the attempts, to show the efficiency of the rejection sampling
        {
            CSV csv;
            SetCSV(csv, GetCSVCols(csv), 0, "Survival");
            WriteSurvival(csv, "white/white", attempts[0], 100);
            WriteSurvival(csv, "white/LDS", attempts[1], 100);
            WriteSurvival(csv, "LDS/white", attempts[2], 100);
            WriteSurvival(csv, "LDS/LDS", attempts[3], 100);
            WriteSurvival(csv, "R2 LDS", attempts[4], 100);
            WriteSurvival(csv, "Sobol LDS", attempts[5], 100);

            SetCSV(csv, GetCSVCols(csv), 0, "Survival Avg");
            WriteSurvival(csv, "white/white", attemptsAvg[0], 100);
            WriteSurvival(csv, "white/LDS", attemptsAvg[1], 100);
            WriteSurvival(csv, "LDS/white", attemptsAvg[2], 100);
            WriteSurvival(csv, "LDS/LDS", attemptsAvg[3], 100);
            WriteSurvival(csv, "R2 LDS", attemptsAvg[4], 100);
            WriteSurvival(csv, "Sobol LDS", attemptsAvg[5], 100);

            SetCSV(csv, GetCSVCols(csv), 0, "Survival Std Dev");
            WriteSurvivalStdDev(csv, "white/white", attemptsAvg[0], attemptsSqAvg[0], 100);
            WriteSurvivalStdDev(csv, "white/LDS", attemptsAvg[1], attemptsSqAvg[1], 100);
            WriteSurvivalStdDev(csv, "LDS/white", attemptsAvg[2], attemptsSqAvg[2], 100);
            WriteSurvivalStdDev(csv, "LDS/LDS", attemptsAvg[3], attemptsSqAvg[3], 100);
            WriteSurvivalStdDev(csv, "R2 LDS", attemptsAvg[4], attemptsSqAvg[4], 100);
            WriteSurvivalStdDev(csv, "Sobol LDS", attemptsAvg[5], attemptsSqAvg[5], 100);

            WriteCSV(csv, "out/uni_cub_survival.csv");
        }

        // make the reports
        for (size_t reportIndex = 0; reportIndex < c_numReportValues; ++reportIndex)
        {
            size_t sampleCount = c_reportValues[reportIndex];

            // open the file
            char buffer[256];
            sprintf_s(buffer, "out/uni_cub_%zu.csv", sampleCount);
            FILE* file = nullptr;
            fopen_s(&file, buffer, "w+t");

            // write the expected histogram
            fprintf(file, "\"Expected\"");
            for (float f : expectedHistogram)
                fprintf(file, ",\"%f\"", f);
            fprintf(file, "\n");

            // write histograms
            WriteHistogram(file, "white/white", histogram[reportIndex][0]);
            WriteHistogram(file, "white/LDS", histogram[reportIndex][1]);
            WriteHistogram(file, "LDS/white", histogram[reportIndex][2]);
            WriteHistogram(file, "LDS/LDS", histogram[reportIndex][3]);
            WriteHistogram(file, "R2 LDS", histogram[reportIndex][4]);
            WriteHistogram(file, "Sobol LDS", histogram[reportIndex][5]);

            // write error
            fprintf(file, "\n\"Error:\"\n");
            WriteHistogramError(file, "white/white", histogram[reportIndex][0], expectedHistogram);
            WriteHistogramError(file, "white/LDS", histogram[reportIndex][1], expectedHistogram);
            WriteHistogramError(file, "LDS/white", histogram[reportIndex][2], expectedHistogram);
            WriteHistogramError(file, "LDS/LDS", histogram[reportIndex][3], expectedHistogram);
            WriteHistogramError(file, "R2 LDS", histogram[reportIndex][4], expectedHistogram);
            WriteHistogramError(file, "Sobol LDS", histogram[reportIndex][5], expectedHistogram);

            // write average error
            fprintf(file, "\n\"Avg Error:\"\n");
            WriteHistogramError(file, "white/white", histogramAvg[reportIndex][0], expectedHistogram);
            WriteHistogramError(file, "white/LDS", histogramAvg[reportIndex][1], expectedHistogram);
            WriteHistogramError(file, "LDS/white", histogramAvg[reportIndex][2], expectedHistogram);
            WriteHistogramError(file, "LDS/LDS", histogramAvg[reportIndex][3], expectedHistogram);
            WriteHistogramError(file, "R2 LDS", histogramAvg[reportIndex][4], expectedHistogram);
            WriteHistogramError(file, "Sobol LDS", histogramAvg[reportIndex][5], expectedHistogram);

            // write std dev
            fprintf(file, "\n\"Std Dev:\"\n");
            WriteHistogramStdDev(file, "white/white", histogramAvg[reportIndex][0], histogramSqAvg[reportIndex][0]);
            WriteHistogramStdDev(file, "white/LDS", histogramAvg[reportIndex][1], histogramSqAvg[reportIndex][1]);
            WriteHistogramStdDev(file, "LDS/white", histogramAvg[reportIndex][2], histogramSqAvg[reportIndex][2]);
            WriteHistogramStdDev(file, "LDS/LDS", histogramAvg[reportIndex][3], histogramSqAvg[reportIndex][3]);
            WriteHistogramStdDev(file, "R2 LDS", histogramAvg[reportIndex][4], histogramSqAvg[reportIndex][4]);
            WriteHistogramStdDev(file, "Sobol LDS", histogramAvg[reportIndex][5], histogramSqAvg[reportIndex][5]);

            fclose(file);
        }
    }

    // uniform to cubic, LDS swap
    {
        // calculate the expected histogram
        std::vector<float> expectedHistogram;
        for (size_t bucketIndex = 0; bucketIndex < c_numHistogramBuckets; ++bucketIndex)
        {
            float p0 = float(bucketIndex) / float(c_numHistogramBuckets);
            float p1 = float(bucketIndex + 1) / float(c_numHistogramBuckets);

            float cdf0 = PDF::Cubic<PDF::UniformWhite, PDF::UniformWhite>::CDF(p0);
            float cdf1 = PDF::Cubic<PDF::UniformWhite, PDF::UniformWhite>::CDF(p1);

            expectedHistogram.push_back(cdf1 - cdf0);
        }

        std::vector<float> histogram[c_numReportValues][6];
        std::vector<float> histogramAvg[c_numReportValues][6];
        std::vector<float> histogramSqAvg[c_numReportValues][6];

        std::vector<size_t> attempts[6];
        std::vector<float> attemptsAvg[6];
        std::vector<float> attemptsSqAvg[6];

        for (size_t testIndex = 0; testIndex < c_numTests; ++testIndex)
        {
            // generate the samples
            std::vector<float> samples[6];
            GenerateSequence<PDF::Cubic<PDF::UniformWhite, PDF::UniformWhite>>(samples[0], c_maxReportValue, 1.0f, attempts[0]);
            GenerateSequence<PDF::Cubic<PDF::UniformWhite, PDF::UniformLDS_Root2>>(samples[1], c_maxReportValue, 1.0f, attempts[1]);
            GenerateSequence<PDF::Cubic<PDF::UniformLDS_GR, PDF::UniformWhite>>(samples[2], c_maxReportValue, 1.0f, attempts[2]);
            GenerateSequence<PDF::Cubic<PDF::UniformLDS_GR, PDF::UniformLDS_Root2>>(samples[3], c_maxReportValue, 1.0f, attempts[3]);
            GenerateSequence<PDF::Cubic<PDF::UniformLDS_R2_x, PDF::UniformLDS_R2_y>>(samples[4], c_maxReportValue, 1.0f, attempts[4]);
            GenerateSequence<PDF::Cubic<PDF::UniformLDS_Sobol_x, PDF::UniformLDS_Sobol_y>>(samples[5], c_maxReportValue, 1.0f, attempts[5]);

            // combine the attempts data
            CombineAttemptsData(attempts[0], attemptsAvg[0], attemptsSqAvg[0], testIndex);
            CombineAttemptsData(attempts[1], attemptsAvg[1], attemptsSqAvg[1], testIndex);
            CombineAttemptsData(attempts[2], attemptsAvg[2], attemptsSqAvg[2], testIndex);
            CombineAttemptsData(attempts[3], attemptsAvg[3], attemptsSqAvg[3], testIndex);
            CombineAttemptsData(attempts[4], attemptsAvg[4], attemptsSqAvg[4], testIndex);
            CombineAttemptsData(attempts[5], attemptsAvg[5], attemptsSqAvg[5], testIndex);

            // calculate data for each report
            for (size_t reportIndex = 0; reportIndex < c_numReportValues; ++reportIndex)
            {
                size_t sampleCount = c_reportValues[reportIndex];

                // calculate the histograms
                CalculateHistogram(samples[0], sampleCount, histogram[reportIndex][0]);
                CalculateHistogram(samples[1], sampleCount, histogram[reportIndex][1]);
                CalculateHistogram(samples[2], sampleCount, histogram[reportIndex][2]);
                CalculateHistogram(samples[3], sampleCount, histogram[reportIndex][3]);
                CalculateHistogram(samples[4], sampleCount, histogram[reportIndex][4]);
                CalculateHistogram(samples[5], sampleCount, histogram[reportIndex][5]);

                // combine the histograms
                HistogramCombine(histogram[reportIndex][0], histogramAvg[reportIndex][0], histogramSqAvg[reportIndex][0], testIndex);
                HistogramCombine(histogram[reportIndex][1], histogramAvg[reportIndex][1], histogramSqAvg[reportIndex][1], testIndex);
                HistogramCombine(histogram[reportIndex][2], histogramAvg[reportIndex][2], histogramSqAvg[reportIndex][2], testIndex);
                HistogramCombine(histogram[reportIndex][3], histogramAvg[reportIndex][3], histogramSqAvg[reportIndex][3], testIndex);
                HistogramCombine(histogram[reportIndex][4], histogramAvg[reportIndex][4], histogramSqAvg[reportIndex][4], testIndex);
                HistogramCombine(histogram[reportIndex][5], histogramAvg[reportIndex][5], histogramSqAvg[reportIndex][5], testIndex);
            }
        }

        // write the attempts, to show the efficiency of the rejection sampling
        {
            CSV csv;
            SetCSV(csv, GetCSVCols(csv), 0, "Survival");
            WriteSurvival(csv, "white/white", attempts[0], 100);
            WriteSurvival(csv, "white/LDS", attempts[1], 100);
            WriteSurvival(csv, "LDS/white", attempts[2], 100);
            WriteSurvival(csv, "LDS/LDS", attempts[3], 100);
            WriteSurvival(csv, "R2 LDS", attempts[4], 100);
            WriteSurvival(csv, "Sobol LDS", attempts[5], 100);

            SetCSV(csv, GetCSVCols(csv), 0, "Survival Avg");
            WriteSurvival(csv, "white/white", attemptsAvg[0], 100);
            WriteSurvival(csv, "white/LDS", attemptsAvg[1], 100);
            WriteSurvival(csv, "LDS/white", attemptsAvg[2], 100);
            WriteSurvival(csv, "LDS/LDS", attemptsAvg[3], 100);
            WriteSurvival(csv, "R2 LDS", attemptsAvg[4], 100);
            WriteSurvival(csv, "Sobol LDS", attemptsAvg[5], 100);

            SetCSV(csv, GetCSVCols(csv), 0, "Survival Std Dev");
            WriteSurvivalStdDev(csv, "white/white", attemptsAvg[0], attemptsSqAvg[0], 100);
            WriteSurvivalStdDev(csv, "white/LDS", attemptsAvg[1], attemptsSqAvg[1], 100);
            WriteSurvivalStdDev(csv, "LDS/white", attemptsAvg[2], attemptsSqAvg[2], 100);
            WriteSurvivalStdDev(csv, "LDS/LDS", attemptsAvg[3], attemptsSqAvg[3], 100);
            WriteSurvivalStdDev(csv, "R2 LDS", attemptsAvg[4], attemptsSqAvg[4], 100);
            WriteSurvivalStdDev(csv, "Sobol LDS", attemptsAvg[5], attemptsSqAvg[5], 100);

            WriteCSV(csv, "out/uni_cub_2_survival.csv");
        }

        // make the reports
        for (size_t reportIndex = 0; reportIndex < c_numReportValues; ++reportIndex)
        {
            size_t sampleCount = c_reportValues[reportIndex];

            // open the file
            char buffer[256];
            sprintf_s(buffer, "out/uni_cub_2_%zu.csv", sampleCount);
            FILE* file = nullptr;
            fopen_s(&file, buffer, "w+t");

            // write the expected histogram
            fprintf(file, "\"Expected\"");
            for (float f : expectedHistogram)
                fprintf(file, ",\"%f\"", f);
            fprintf(file, "\n");

            // write histograms
            WriteHistogram(file, "white/white", histogram[reportIndex][0]);
            WriteHistogram(file, "white/LDS", histogram[reportIndex][1]);
            WriteHistogram(file, "LDS/white", histogram[reportIndex][2]);
            WriteHistogram(file, "LDS/LDS", histogram[reportIndex][3]);
            WriteHistogram(file, "R2 LDS", histogram[reportIndex][4]);
            WriteHistogram(file, "Sobol LDS", histogram[reportIndex][5]);

            // write error
            fprintf(file, "\n\"Error:\"\n");
            WriteHistogramError(file, "white/white", histogram[reportIndex][0], expectedHistogram);
            WriteHistogramError(file, "white/LDS", histogram[reportIndex][1], expectedHistogram);
            WriteHistogramError(file, "LDS/white", histogram[reportIndex][2], expectedHistogram);
            WriteHistogramError(file, "LDS/LDS", histogram[reportIndex][3], expectedHistogram);
            WriteHistogramError(file, "R2 LDS", histogram[reportIndex][4], expectedHistogram);
            WriteHistogramError(file, "Sobol LDS", histogram[reportIndex][5], expectedHistogram);

            // write average error
            fprintf(file, "\n\"Avg Error:\"\n");
            WriteHistogramError(file, "white/white", histogramAvg[reportIndex][0], expectedHistogram);
            WriteHistogramError(file, "white/LDS", histogramAvg[reportIndex][1], expectedHistogram);
            WriteHistogramError(file, "LDS/white", histogramAvg[reportIndex][2], expectedHistogram);
            WriteHistogramError(file, "LDS/LDS", histogramAvg[reportIndex][3], expectedHistogram);
            WriteHistogramError(file, "R2 LDS", histogramAvg[reportIndex][4], expectedHistogram);
            WriteHistogramError(file, "Sobol LDS", histogramAvg[reportIndex][5], expectedHistogram);

            // write std dev
            fprintf(file, "\n\"Std Dev:\"\n");
            WriteHistogramStdDev(file, "white/white", histogramAvg[reportIndex][0], histogramSqAvg[reportIndex][0]);
            WriteHistogramStdDev(file, "white/LDS", histogramAvg[reportIndex][1], histogramSqAvg[reportIndex][1]);
            WriteHistogramStdDev(file, "LDS/white", histogramAvg[reportIndex][2], histogramSqAvg[reportIndex][2]);
            WriteHistogramStdDev(file, "LDS/LDS", histogramAvg[reportIndex][3], histogramSqAvg[reportIndex][3]);
            WriteHistogramStdDev(file, "R2 LDS", histogramAvg[reportIndex][4], histogramSqAvg[reportIndex][4]);
            WriteHistogramStdDev(file, "Sobol LDS", histogramAvg[reportIndex][5], histogramSqAvg[reportIndex][5]);

            fclose(file);
        }
    }

    return 0;
}
