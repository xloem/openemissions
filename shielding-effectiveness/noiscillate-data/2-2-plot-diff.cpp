// 2-2-plot-diff.cpp
//
// Produces charts of isolated emitter power (given a recording with and without emitter)
// and shielding effectiveness betweeen two environments (given four recordings)

#include <memory>
#include <string>

#include <TCanvas.h>

#include "chartprofile.hpp"
#include "processornoiseprofile.hpp"
#include "rootapplication.hpp"

using Scalar = double;

class T22PlotDiff : public RootApplication
{
public:
	using NoiseP = ProcessorNoiseProfile<Scalar, StatsAccumulatorHistogram<Scalar>, STATS_VARIANCE>;

	T22PlotDiff(Int_t * argc, char ** argv)
	: RootApplication(
		"2-2-plot-diff",
		"Calculates changes between emitters and environments to estimate either\n"
		"the power spectrum of one emitter, or the shielding effectiveness between\n"
		"two otherwise identical environments.",
		argc, argv),
	  _sePad(nullptr)
	{
		regArg({{"-u","--unshielded"}, {"bg.noisep"}, {"fg.noisep"}}, {"profiles to compare", "fg.noisep must have 1 more emitter than bg.noisep"},
			[this](std::string, std::string bgfn, std::string fgfn)
			{
				_unshielded.reset(new ChartedEmitterProfile(bgfn, fgfn));
			}, 1 // required
		);
		regArg({{"-s","--shielded"}, {"bg.noisep"}, {"fg.noisep"}}, {"additional comparison for shielding effectiveness", "fg.noisep must have 1 more emitter than bg.noisep"},
			[this](std::string, std::string bgfn, std::string fgfn)
			{
				_shielded.reset(new ChartedEmitterProfile(bgfn, fgfn));
			}, 0 // not required
		);
	}

	virtual void GetOptions(Int_t * argc, char **argv) override
	{
		RootApplication::GetOptions(argc, argv);
		if (_shielded)
		{
			if (_shielded->emitter() != _unshielded->emitter())
			{
				std::cerr << "Emitter of interest differs across recording pairs." << std::endl;
				Terminate(-1);
			}
		}
		if (quiet < 1)
		{
			_canvas.reset(new TCanvas());
			
			if (_shielded)
			{
				_canvas->cd();
				_canvas->Divide(1,2);
				_canvas->GetPad(1)->Divide(2,2);
				_unshielded->setPad(_canvas->GetPad(1)->GetPad(1));
				_shielded->setPad(_canvas->GetPad(1)->GetPad(2));
				_sePad = _canvas->GetPad(1);
			}
			else
			{
				_unshielded->setPad(&*_canvas);
			}
		}
	}

	virtual void Run(Bool_t retrn) override
	{
		_unshielded->paint();
		if (_shielded)
		{
			_shielded->paint();
		}
		tick();
		if (_shielded)
		{
			_seP.set(_shielded->profile(), _unshielded->profile());
			_seChart.reset(new RootChartProfile<BriefDiffProfile>("Shielding Effectiveness"));
			_seChart->pad() = _sePad;
			_seChart->prepPoints(_seP, quiet > 0);
			_seChart->paint(_seP);
			_seChart->finalize();
			_seChart->paint(_seP);
			tick();
		}
	}

private:
	class BriefDiffProfile
	{
	public:
		using Scalar = ::Scalar;
		decltype(auto) begin() { return _freqs.cbegin(); }
		decltype(auto) end() { return _freqs.cend(); }
		Scalar bestMag(Scalar freq) { return _dists.at(freq).mean(); }
		Scalar bestMagVariance(Scalar freq) { return _dists.at(freq).variance(); }
		size_t size() { return _freqs.size(); }

		template <typename Profile>
		void set(Profile & bg, Profile & fg)
		{
			_dists.clear();
			_freqs.clear();
			for (auto freq : bg)
			{
				try
				{
					auto mag = fg.bestMag(freq) - bg.bestMag(freq);
					auto var = fg.bestMagVariance(freq) - bg.bestMagVariance(freq);
					_freqs.push_back(freq);
					// funny syntax provides for true map emplacing
					_dists.emplace(
						std::piecewise_construct,
						std::forward_as_tuple(freq),
						std::forward_as_tuple(mag, var)
					);
				}
				catch (std::out_of_range)
				{ }
			}
		}
	private:
		std::vector<Scalar> _freqs;
		std::map<Scalar, StatsDistributionBrief<Scalar>> _dists;
	};
	class ChartedEmitterProfile
	{
	public:
		ChartedEmitterProfile(std::string bgfn, std::string fgfn, TVirtualPad * pad = nullptr)
		: _bgChart("No Emitter"),
		  _fgChart("With Emitter"),
		  _emChart("Emitter"),
		  _bgP(bgfn, {0}),
		  _fgP(fgfn, {0}, _bgP.environment(), {}, _bgP.srcDescs())
		{
    			// find emitter from difference
			std::set<int64_t> bgE{_bgP.emitters().begin(), _bgP.emitters().end()};
			std::set<int64_t> fgE{_fgP.emitters().begin(), _fgP.emitters().end()};

			std::cerr << "Emitters in " << bgfn << ":";
			for (auto emitter : bgE)
			{
			  std::cerr << " " << emitter;
			}
			std::cerr << std::endl;
			std::cerr << "Emitters in " << fgfn << ":";
			for (auto emitter : fgE)
			{
			  std::cerr << " " << emitter;
			}
			std::cerr << std::endl;

			// symmetric difference of two sets
			std::vector<int64_t> emitterVec(bgE.size()+fgE.size());
			auto tail = std::set_symmetric_difference(bgE.begin(), bgE.end(), fgE.begin(), fgE.end(), emitterVec.begin());
			emitterVec.resize(tail - emitterVec.begin());

			if (emitterVec.size() != 1)
			{
			  throw std::invalid_argument("Profiles must differ by exactly 1 active emitter.");
			}

			_emitter = emitterVec[0];
			if (bgE.count(_emitter))
			{
			  throw std::invalid_argument("Profiles provided in reverse order.");
			}

			std::cerr << "Producing a spectrum of emitter " << _emitter << " by calculating " << fgfn << " - " << bgfn << std::endl;

			_emP.set(_bgP, _fgP);
		}
		void setPad(TVirtualPad * pad)
		{
			if (pad != nullptr)
			{
				pad->cd();
				pad->Divide(1,2);
				pad->GetPad(1)->Divide(2,1);
				_bgChart.pad() = pad->GetPad(1)->GetPad(1);
				_fgChart.pad() = pad->GetPad(1)->GetPad(2);
				_emChart.pad() = pad->GetPad(2);
			}
			_bgChart.prepPoints(_bgP, !pad);
			_fgChart.prepPoints(_fgP, !pad);
			_emChart.prepPoints(_emP, !pad);
		}
		void paint()
		{
			_bgChart.paint(_bgP);
			_fgChart.paint(_fgP);
			_bgChart.finalize();
			_fgChart.finalize();
			_bgChart.paint(_bgP);
			_fgChart.paint(_fgP);

			_emChart.paint(_emP);
			_emChart.finalize();
			_emChart.paint(_emP);
		}
		int64_t emitter() const { return _emitter; }
		BriefDiffProfile & profile() { return _emP; }
	private:
		RootChartProfile<NoiseP> _bgChart;
		RootChartProfile<NoiseP> _fgChart;
		RootChartProfile<BriefDiffProfile> _emChart;
		int64_t _emitter;

		NoiseP _bgP, _fgP;
		BriefDiffProfile _emP;
	};

	std::unique_ptr<ChartedEmitterProfile> _unshielded, _shielded;

	std::unique_ptr<TCanvas> _canvas;
	TVirtualPad * _sePad;

	BriefDiffProfile _seP;
	std::unique_ptr<RootChartProfile<BriefDiffProfile>> _seChart;
};

int main(int argc, char ** argv)
{
  T22PlotDiff app(&argc, argv);
  app.GetOptions(&argc, argv);
  app.Run(false);
  return 0;
}
