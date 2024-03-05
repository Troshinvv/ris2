#ifndef BUNCH_H
#define BUNCH_H

#include "wrap.h"
#include <Axis.hpp>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

class Palette{
public:
  Palette() = default;
  ~Palette() = default;
  
  Palette& SetPalette( std::vector<Style> palette ){ styles_ = std::move(palette); return *this; }
  const std::vector<Style>& GetPalette(){ return styles_; }
  template< class T >
  void PaintObjects( std::vector< Wrap<T> >& objects ) const { 
    std::for_each( objects.begin(), objects.end(), 
    [this, i=0]( auto& obj ) mutable {
      obj.SetStyle( styles_.at(i) );
      ++i;
    } ); 
  }
  void PaintObject( size_t style_idx, TGraphErrors* graph ) { 
    styles_.at(style_idx).operator()( graph );
  }
  TLegend* MakeLegend( const std::vector<std::string>& captions, std::vector<double> position = {} ){
    int i=0;
    TLegend* leg{};
    if( position.size() == 4 )
      leg = new TLegend{ position.at(0), position.at(1), position.at(2), position.at(3) };
    else{
      leg = new TLegend{};
    }
    for( auto& titile : captions ){
      auto graph = new TGraph(0);
      graph->SetLineColor( styles_.at(i).color_ );
      graph->SetMarkerColor( styles_.at(i).color_ );
      if( styles_.at(i).marker_ >= 0 ){
        graph->SetMarkerStyle( styles_.at(i).marker_ );
        leg->AddEntry( graph, titile.c_str(), "P" );
      }
      if( styles_.at(i).marker_ < 0 ){
        graph->SetLineStyle( styles_.at(i).marker_ );
        leg->AddEntry( graph, titile.c_str(), "L" );
      }
      ++i;
    }
    return leg;
  }
protected:
  std::vector<Style> styles_{};
};

template<typename T>
class Bunch {
public:
  Bunch<T>() = default;
  Bunch<T>( const Bunch<T>& other) = default;
  Bunch<T>( Bunch<T>&& other) = default;
  Bunch<T>& operator=( const Bunch<T>& other) = default;
  Bunch<T>& operator=( Bunch<T>&& ) = default;

  ~Bunch<T>() = default;

  template<typename... Args>
  Bunch& AddToBunch( std::string title, Args... args ){ 
    bunch_.emplace_back( title, args... ); 
    return *this;
  }
  Palette& GetPalette() { return palette_; }
  template<typename Func>
  Bunch& Perform( Func function ){ std::for_each( bunch_.begin(), bunch_.end(), function ); return *this; }
  std::vector<Wrap<T>>& operator*(){
    palette_.PaintObjects(bunch_);
    return bunch_; 
  }
  template<typename Func>
  Bunch& operator()( Func function ){ std::for_each( bunch_.begin(), bunch_.end(), function ); return *this; }
  Wrap<T>& operator[]( size_t idx ){ return bunch_.at(idx); }
  TLegend* MakeLegend(std::vector<double> position = {}){
    std::vector<std::string> bunch_titles_;
    std::for_each( bunch_.begin(), bunch_.end(), [&bunch_titles_]( const Wrap<T>& obj ) mutable { 
      bunch_titles_.emplace_back( obj.GetTitle() ); 
    } );
    auto leg = palette_.MakeLegend( bunch_titles_, position );
    return leg;
  }

private:
  Palette palette_{};
  std::vector< Wrap<T> > bunch_{};
};

class DoubleDifferential {
public:
  DoubleDifferential() = delete;
  template<typename ...Args>
  DoubleDifferential(std::string title, Args... args) : base_correlation_{ title, args... }{ }
  DoubleDifferential( const DoubleDifferential& ) = default;
  DoubleDifferential& operator=( const DoubleDifferential& ) = default;
  DoubleDifferential( DoubleDifferential&& ) = default;
  DoubleDifferential& operator=( DoubleDifferential&& ) = default;
  ~DoubleDifferential() = default;
  template<typename Func>
  DoubleDifferential& Perform( const Func& function ){ 
    function( base_correlation_ ); 
    return *this; 
  }
  template<typename Func>
  DoubleDifferential& operator()( const Func& function ){ 
    function( base_correlation_ ); 
    return *this; 
  }
  Palette& GetPalette(){ return palette_; }
  DoubleDifferential& SetSliceAxis(Qn::AxisD slie_axis){ slice_axis_ = std::move(slie_axis); return *this; }
  DoubleDifferential& SetProjectionAxis(Qn::AxisD projection_axis){ projection_axis_ = std::move(projection_axis); return *this; }
  DoubleDifferential& SetPalette( const std::vector<Style>& styles ){
    palette_.SetPalette(styles);
    return *this;
  }
  TLegend* MakeLegend(const std::string& slice_var, const std::string& units="", std::vector<double> position = {} ){
    auto bin_edges = slice_axis_.GetBinEdges();
    std::vector<std::string> captions{};
    for( int i = 0; i < bin_edges.size() - 1; ++i ){
      auto stream = std::ostringstream{};
      stream << std::setprecision(2) << slice_var << " " << bin_edges.at(i) << "-" << bin_edges.at(i+1) << " " << units;
      captions.emplace_back( stream.str() );
    }
    return palette_.MakeLegend( captions, position );
  }
  std::vector<Wrap<Correlation>>& operator*(){
    if( result_correlations_.empty() )
      Slice();
    return result_correlations_;
  }
private:
  void Slice(){
    auto n_projections = slice_axis_.GetNBins();
    auto slice_var = slice_axis_.Name();
    auto bin_edges = slice_axis_.GetBinEdges();
    // auto projection_var = projection_axis_.Name();
    base_correlation_
      .Rebin( std::vector<Qn::AxisD>{slice_axis_, projection_axis_} )
      .Project( std::vector<Qn::AxisD>{slice_axis_, projection_axis_} );
    result_correlations_.reserve(n_projections);
    for( int i=0; i<n_projections; ++i ){
      result_correlations_.emplace_back( base_correlation_ );
      
      result_correlations_.back().Rebin( std::vector<Qn::AxisD>{{slice_var, 1, bin_edges.at(i), bin_edges.at(i+1)} } )
        .Project( std::vector<Qn::AxisD>{{projection_axis_}} );
      palette_.PaintObjects( result_correlations_ );
    }
  }
  Palette palette_;
  Wrap<Correlation> base_correlation_;
  std::vector< Wrap<Correlation> > result_correlations_;
  Qn::AxisD slice_axis_{};
  Qn::AxisD projection_axis_{};
};

template<typename T>
class RatioBuilder {
public:
  template<typename... Args>
  RatioBuilder( std::string title, Args... args ) {
    results_.emplace_back( args... );
    titles_.emplace_back( std::move(title) );
  }
  ~RatioBuilder() = default;
  template<typename ...Args>
  RatioBuilder& AddToBunch( std::string title, Args... args ){ 
    titles_.emplace_back( title );
    results_.emplace_back( args...); 
    return *this; 
  }
  std::vector<Result<T>>& GetResults(){
    return results_;
  }
  std::vector<Result<T>>& GetRatios(){
    BuildRatios();
    return ratios_;
  }
  template<typename Func>
  RatioBuilder<T>& Perform( const Func& function ){ 
    std::for_each( results_.begin(), results_.end(), function ); 
    return *this; 
  }
  template<typename Func>
  RatioBuilder<T>& operator()( const Func& function ){ 
    std::for_each( results_.begin(), results_.end(), function ); 
    return *this; 
  }
  RatioBuilder<T>& SetPalette( std::vector<Style> corr_styles ){
    palette_.SetPalette(corr_styles);
    return *this;
  }
  Result<T>& GetReference(){ return results_.front(); }
  TLegend* MakeLegend( std::vector<double> position = {} ){
    auto leg = palette_.MakeLegend(titles_, position);
    return leg;
  }
  size_t Size() const { return results_.size(); }
  std::vector<Wrap<TGraphErrors>>& GetResultWraps(){ UpdatePoints(); return result_points_; }
  std::vector<Wrap<TGraphErrors>>& GetRatioWraps(){ UpdatePoints(); return ratio_points_; }
private:
  void BuildRatios(){
    if( ratios_.empty() ){
      int i=0;
      const auto& reference = results_.front();
      ratios_.reserve( results_.size() );
      std::for_each( results_.begin(), results_.end(), [&reference, this, i=0](auto& obj) mutable {
        ratios_.emplace_back( obj.Divide( reference )  );
      } );
    }
  }
  void UpdatePoints(){
    BuildRatios();
    if( ! result_points_.empty() ){
      return;
    }
    std::for_each( results_.begin(), results_.end(), 
    [this]( auto obj  ) {
      result_points_.push_back( Wrap<TGraphErrors>(obj.GetPoints()) );
    });
    std::for_each( ratios_.begin(), ratios_.end(), 
    [this]( auto obj  ) {
      ratio_points_.push_back( Wrap<TGraphErrors>(obj.GetPoints()) );
    });
    palette_.PaintObjects(result_points_);
    palette_.PaintObjects(ratio_points_);
  };
  std::vector<std::string>  titles_{};
  Palette palette_;
  std::vector<Result<T>> results_{};
  std::vector<Result<T>> ratios_{};
  std::vector<Wrap<TGraphErrors>> result_points_{};
  std::vector<Wrap<TGraphErrors>> ratio_points_{};
};

}

#endif // BUNCH_H
