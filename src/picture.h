#ifndef PICTURE_H
#define PICTURE_H

#include <TPad.h>
#include <TVirtualPad.h>
#include <algorithm>
#include <cstddef>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TLatex.h>
#include <TF1.h>
#include <TLegend.h>

#include "wrap.h"
#include "bunch.h"

struct Axis{
  Axis() = default;
  Axis& SetTitle( std::string title ){ title_ = std::move(title); return *this; }
  Axis& SetHi( double hi ){ hi_ = hi; return *this; }
  Axis& SetLo( double lo ){ lo_ = lo; return *this; }
  std::string title_{};
  double hi_{};
  double lo_{};
};

struct Text{
  Text() = default;
  Text& SetSize(double size){ size_ = size; return *this; }
  Text& SetText(std::string text){ text_ = std::move(text); return *this; }
  Text& SetPosition(std::array<double, 2> pos){ position_ = std::move(pos); return *this; }
  double size_;
  std::string text_;
  std::array<double, 2> position_;
};

class Picture{
public:
  Picture() : graph_stack_{ std::make_unique<TMultiGraph>() } {}
  Picture(const std::vector<double>& aspect_ratio ) : graph_stack_{ std::make_unique<TMultiGraph>() } {
    auto rd = std::random_device{};
    auto re = std::mt19937{rd()};
    auto distr = std::uniform_int_distribution<size_t>{ 0, std::numeric_limits<size_t>::max() };
    auto name = std::string{"canv_" }.append( std::to_string(distr(re)) );
    pad_ = new TPad( name.c_str(), "", 
                                    aspect_ratio.at(0), 
                                    aspect_ratio.at(1), 
                                    aspect_ratio.at(2), 
                                    aspect_ratio.at(3) );
  }
  Picture( Picture&& ) = default;
  Picture( const Picture& ) = delete;
  template<typename T>
  Picture& AddToPlot( Wrap<T>& wrap ){ 
    auto style = wrap.GetStyle();
    if( style.marker_ >= 0 )
      graph_stack_->Add( wrap.ReleasePoints(), "P" ); 
    if( style.marker_ < 0 )
      graph_stack_->Add( wrap.ReleasePoints(), "L" ); 
    return *this; 
  }
  Picture& SetResolution( std::vector<double> aspect_ratio ){
    auto rd = std::random_device{};
    auto re = std::mt19937{rd()};
    auto distr = std::uniform_int_distribution<size_t>{ 0, std::numeric_limits<size_t>::max() };
    auto name = std::string{"pad_" }.append( std::to_string(distr(re)) );
    pad_ = new TPad( name.c_str(), "", 
                              aspect_ratio.at(0), 
                              aspect_ratio.at(1), 
                              aspect_ratio.at(2), 
                              aspect_ratio.at(3) );
    return *this;
  }
  Picture& SetXAxis( Axis axis ){ x_axis_ = axis; return *this; }
  Picture& SetYAxis( Axis axis ){ y_axis_ = axis; return *this; }
  template<typename T>
  Picture& AddToPlot( std::vector<Wrap<T>>& bunch ){
    std::for_each( bunch.begin(), bunch.end(), [this]( Wrap<T>& obj ){ 
      if( obj.GetStyle().marker_ >= 0 )
        graph_stack_->Add( obj.ReleasePoints(), "P" ); 
      if( obj.GetStyle().marker_ < 0 )
        graph_stack_->Add( obj.ReleasePoints(), "L" ); 
    } );    
    return *this;
  }
  Picture& AddText( Text text ) {
    texts_.emplace_back(std::move(text));
    return *this;
  }
  Picture& AddLegend( TLegend* leg ){
    legends_.emplace_back( leg );
    return *this;
  }
  TPad* Print(  ){
    pad_->cd();
    auto axis_titles = std::string{";"}.append( x_axis_.title_ ).append(";").append(y_axis_.title_);
    graph_stack_->SetTitle(axis_titles.c_str());
    graph_stack_->Draw("APL");
    if( x_axis_.hi_ > x_axis_.lo_ ){
      graph_stack_->GetHistogram()->GetXaxis()->SetRangeUser(x_axis_.lo_, x_axis_.hi_);
      graph_stack_->GetHistogram()->Draw("APL");
    }
    if( y_axis_.hi_ > y_axis_.lo_ ){
      graph_stack_->GetHistogram()->GetYaxis()->SetRangeUser(y_axis_.lo_, y_axis_.hi_);
      graph_stack_->Draw("APL");
    }
    vec_tlatexs_.reserve( texts_.size() );
    for( auto text : texts_ ){
      vec_tlatexs_.emplace_back( new TLatex( text.position_.at(0), text.position_.at(1), text.text_.c_str() ) );
      vec_tlatexs_.back()->SetNDC();
      vec_tlatexs_.back()->SetTextSize( text.size_ );
      vec_tlatexs_.back()->Draw("same");
    }
    for( auto& leg : legends_ ){
      leg->Draw( "same" );
    }
    zero_line_.reset(new TF1( "zero", "0", x_axis_.lo_, x_axis_.hi_ ));
    zero_line_->Draw("same");
    return pad_;
  }
private:
  Axis x_axis_;
  Axis y_axis_;
  TPad* pad_{};
  std::unique_ptr<TMultiGraph> graph_stack_{};
  std::vector<Text> texts_{};
  std::vector< std::unique_ptr<TLatex> > vec_tlatexs_{};
  std::unique_ptr<TF1> zero_line_{};
  std::vector<std::unique_ptr<TLegend>> legends_;
};


class Plot{
public:
  Plot( std::vector<size_t> resolution ){
    auto rd = std::random_device{};
    auto re = std::mt19937{rd()};
    auto distr = std::uniform_int_distribution<size_t>{ 0, std::numeric_limits<size_t>::max() };
    auto name = std::string{"canv_" }.append( std::to_string(distr(re)) );
    canvas_.reset( new TCanvas(name.c_str(), "", resolution.at(0), resolution.at(1)) );
  }
  ~Plot() = default;
  template<typename... Args>
  Picture& AddSubPlot( Args... args ){ 
    plots_.emplace_back(args...); return plots_.back(); 
  }
  Picture& LastSubPlot(){
    if( plots_.empty() ){
      plots_.emplace_back( std::vector<double>{0., 0., 1., 1.} );
    }
    return plots_.back(); 
  }
  Picture& GetSubPlot( size_t idx ){ return plots_.at(idx); }
  template<typename T>
  Plot& AddRatioPlot( RatioBuilder<T>& ratio, 
    std::vector<double> result_plot = { 0.0, 0.35, 1.0, 1.0 },
    std::vector<double> ratio_plot = { 0.0, 0.0, 1.0, .35 }){
      plots_.emplace_back( result_plot );
      plots_.back().AddToPlot( ratio.GetResults() );
      plots_.back().AddToPlot( ratio.GetReference() );
      plots_.emplace_back( ratio_plot );
      plots_.back().AddToPlot( ratio.GetRatios() );
      return *this;
  }
  Plot& Print( const std::string& save_name ){
    canvas_->cd();
    std::vector<TPad*> pads;
    pads.reserve(plots_.size());
    for( auto& plot : plots_ ){
      pads.push_back(plot.Print());
    }
    canvas_->cd();
    for( auto pad : pads ){
      pad->Draw();
    }
    canvas_->Draw();
    canvas_->Print( save_name.c_str() );
    return *this;
  }
protected:
  std::vector<Picture> plots_;
  std::unique_ptr<TCanvas> canvas_{};
};

class RatioPlot : public Plot{
  
};

#endif // PICTURE_H