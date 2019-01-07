#pragma once

#include <initializer_list>
#include <iostream>

#include <TApplication.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TSysEvtHandler.h>
#include <TSystem.h>

class RootApplication : public TApplication
{
public:
  RootApplication(char const * name, std::string helpDescr, Int_t * argc, char ** argv)
  : TApplication(name, argc, argv, 0, -1),
    terminate(false),
    quiet(0),
    _helpDescr(helpDescr),
    _sigQuit(*this, kSigQuit),
    _sigInterrupt(*this, kSigInterrupt),
    _sigTermination(*this, kSigTermination)
  {
    graph.SetTitle(name);
    graph.GetXaxis()->SetTitle("Frequency (Hz)");
    graph.GetYaxis()->SetTitle("Power (dB)");
    graph.SetEditable(kFALSE);
    graph.SetHighlight(kTRUE);

    _sigQuit.Add();
    _sigInterrupt.Add();
    _sigTermination.Add();
  }

protected:
  void tick()
  {
    gSystem->ProcessEvents();
  }

  void regArg(std::initializer_list<std::initializer_list<std::string>> matches, std::vector<std::string> && descr, std::function<void(std::vector<std::string> &)> handler, size_t subargs, size_t requiredCt)
  {
    if (matches.size() != subargs) throw std::logic_error("doc mismatches implementation");
    auto argIdx = _argDescrs.size();
    _argDescrs.emplace_back();
    Arg & arg = _argDescrs.back();
    for (auto & match : matches)
    {
      arg.matches.emplace_back(match);
    }
    arg.handler = handler;
    arg.args = subargs;
    arg.required = requiredCt;
    arg.requiredLeft = requiredCt;
    arg.descr = std::move(descr);
    for (auto & match : arg.matches[0])
    {
      _argMap[match] = argIdx;
    }
  }

  void regArg(std::initializer_list<std::initializer_list<std::string>> matches, std::vector<std::string> descr, std::function<void(std::string)> handler, size_t requiredCt = 0)
  {
    regArg(matches, std::move(descr), [handler](std::vector<std::string> & v) { handler(v[0]); }, 1, requiredCt);
  }

  void regArg(std::initializer_list<std::initializer_list<std::string>> matches, std::vector<std::string> descr, std::function<void(std::string, std::string)> handler, size_t requiredCt = 0)
  {
    regArg(matches, std::move(descr), [handler](std::vector<std::string> & v) { handler(v[0], v[1]); }, 2, requiredCt);
  }

  void regArg(std::initializer_list<std::initializer_list<std::string>> matches, std::vector<std::string> descr, std::function<void(std::string, std::string, std::string)> handler, size_t requiredCt = 0)
  {
    regArg(matches, std::move(descr), [handler](std::vector<std::string> & v) { handler(v[0], v[1], v[2]); }, 3, requiredCt);
  }

  virtual void GetOptions(Int_t *argc, char **argv) override
  {
    regArg({{"-h","-?","--help"}}, {"display this message"},
      [this](std::string)
      {
        PrintHelp();
        Terminate(0);
      }
    );
    regArg({{"-q"}}, {"do not display gfx or show progress"},
      [this](std::string)
      {
        ++ quiet;
      }
    );

    std::vector<std::string> subargs;
    size_t otherCt = 0;
    for (int i = 1; i < *argc; ++ i)
    {
      try
      {
        try
        {
          auto & handler = _argDescrs[_argMap.at(argv[i])];
          subargs.resize(handler.args);
          for (int j = 0, i2 = i; j < handler.args; ++ j, ++i2)
          {
            if (i2 >= *argc)
            {
              throw std::invalid_argument(handler.matches[0][0] + " missing " + handler.matches[j+1][0]);
            }
            subargs[j] = argv[i2];
            argv[i2] = nullptr;
            i = i2;
          }
          handler.handler(subargs);
          if (handler.requiredLeft)
          {
            -- handler.requiredLeft;
          }
          goto continue_outer;
        }
        catch (std::out_of_range)
        {
          size_t otherMatch = 0;
          Arg * lastPossible = 0;
          for (auto & argDescr : _argDescrs)
          {
            if (argDescr.matches[0].size() == 0)
            {
              if (otherMatch == otherCt)
              {
                subargs.resize(2);
                subargs[0] = "";
                subargs[1] = argv[i];
                argDescr.handler(subargs);
                if (argDescr.requiredLeft)
                {
                  -- argDescr.requiredLeft;
                }
                if (argDescr.requiredLeft == 0)
                {
                  ++ otherCt;
                  lastPossible = &argDescr;
                }
                goto continue_outer;
              }
              else
              {
                ++ otherMatch;
              }
            }
          }
          if (lastPossible)
          {
            lastPossible->handler(subargs);
          }
          else
          {
            throw std::invalid_argument(argv[i]);
          }
        }
      }
      catch(std::invalid_argument e)
      {
        std::cerr << "Invalid argument " << e.what() << " in";
      }
      catch(std::out_of_range e)
      {
        std::cerr << "Value out of range " << e.what() << " in";
      }
      for (auto & subarg : subargs)
      {
        std::cerr << " " << subarg;
      }
      PrintHelp();
      argv[i] = nullptr;
      Terminate(-1);
continue_outer:
      ;
    }
    bool success = true;
    for (auto & argDescr : _argDescrs)
    {
      if (argDescr.requiredLeft > 0)
      {
        std::cerr << "Missing required argument:";
        for (auto & match : argDescr.matches)
        {
          if (match.size() > 0)
          {
            std::cerr << " " << match[0];
          }
        }
        std::cerr << std::endl;
        success = false;
      }
    }
    if (success == false)
    {
      PrintHelp();
      Terminate(-1);
    }
  }

  virtual void PrintHelp()
  {
    std::cout << std::endl;
    std::cout << _helpDescr << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: " << ApplicationName() << " [args ...]";
    for (auto & argDescr : _argDescrs)
    {
      auto count = argDescr.required;
      if (count == 0 && argDescr.matches[0].size() == 0)
      {
        count = 1;
        std::cout << " [";
      }
      for (size_t i = 0; i < count; ++ i)
      {
        for (auto & match : argDescr.matches)
        {
          if (match.size())
          {
            if (!(i == 0 && count != argDescr.required))
            {
              // open bracket doesn't precede us (todo make clearer with e.g. flagSpacePrefix set aboev)
              std::cout << " ";
            }
            std::cout << match[0];
          }
        }
      }
      if (count == 1 && argDescr.required == 0)
      {
        std::cout << "]";
      }
    }
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "Arguments:" << std::endl;
    size_t longest = 0;
    for (auto & argDescr : _argDescrs)
    {
      size_t size = 0;
      for (auto & grp : argDescr.matches)
      {
        for (auto & match : grp)
        {
          size += match.size() + 1;
        }
      }
      if (size > longest)
      {
        longest = size;
      }
    }
    for (auto & argDescr : _argDescrs)
    {
      size_t size = 0;
      for (auto & grp : argDescr.matches)
      {
        for (auto & match : grp)
        {
          size += match.size() + 1;
        }
      }
      for (size_t i = 0; i < longest - size; ++ i)
      {
        std::cout << " ";
      }
      for (size_t i = 0; i < argDescr.matches.size(); ++ i)
      {
        for (size_t j = 0; j < argDescr.matches[i].size(); ++ j)
        {
          if (j > 0)
          {
            std::cout << (i > 0 ? "|" : ",");
          }
          else
          {
            std::cout << " ";
          };
          std::cout << argDescr.matches[i][j];
        }
      }
      std::cout << ": " << argDescr.descr[0] << std::endl;
      for (size_t i = 1; i < argDescr.descr.size(); ++ i)
      {
        for (size_t j = 0; j < longest; ++ j)
        {
          std::cout << " ";
        }
        std::cout << "  " << argDescr.descr[i] << std::endl;
      }
    }
    std::cout << std::endl;
  }


  bool terminate;
  int quiet;
  std::unique_ptr<TCanvas> canvas;
  TGraphAsymmErrors graph;

private:
  virtual void HandleException(Int_t sig = -1) override
  {
    std::cerr << "-exception-" << std::endl;
    terminate = true;
  }

  class THaltSignalHandler : public TSignalHandler
  {
  public:
    THaltSignalHandler(RootApplication & target, ESignals sig)
    : TSignalHandler(sig),
      target(target)
    { }

    virtual Bool_t Notify() override
    {
      target.HandleException(GetSignal());
      return true;
    }

  private:
    RootApplication & target;
  };

  struct Arg
  {
    std::vector<std::vector<std::string>> matches;
    std::function<void(std::vector<std::string> &)> handler;
    size_t args;
    size_t required;
    size_t requiredLeft;
    std::vector<std::string> descr;
  };

  std::string _helpDescr;
  std::vector<Arg> _argDescrs;
  std::map<std::string, size_t> _argMap;

  THaltSignalHandler _sigQuit, _sigInterrupt, _sigTermination;

};
