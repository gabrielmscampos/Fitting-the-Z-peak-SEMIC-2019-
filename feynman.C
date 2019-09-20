{
//=========Macro generated from canvas: c1/A canvas
//=========  (Thu Jun 12 12:21:17 2003) by ROOT version3.03/09
   TCanvas *c1 = new TCanvas("c1", "A canvas",0,0,1199,789);
   c1->Range(0,0,140,60);
   c1->SetBorderSize(2);
   c1->SetFrameFillColor(0);

   // TLine *line = new TLine(10,10,30,30);
   // line->SetLineWidth(3);
   // line->Draw();
   // line = new TLine(10,50,30,30);
   // line->SetLineWidth(3);
   // line->Draw();

   TArrow *ar1 = new TArrow(10,10,30,30,0.03, "-|>-");
   ar1->SetLineWidth(3);
   ar1->SetFillColor(1);
   ar1->Draw();

   TArrow *ar2 = new TArrow(10,50,30,30,0.03, "-<|-");
   ar2->SetLineWidth(3);
   ar2->SetFillColor(1);
   ar2->Draw();

   TLatex *tex = new TLatex(7,6,"#font[12]{#bar{q}}");
   tex->SetTextAlign(22);
   tex->SetTextSize(0.1);
   tex->SetLineWidth(2);
   tex->Draw();

      tex = new TLatex(7,55,"#font[12]{q}");
   tex->SetTextAlign(22);
   tex->SetTextSize(0.1);
   tex->SetLineWidth(2);
   tex->Draw();

   //    tex = new TLatex(42.5,37.7,"#gamma");
   // tex->SetTextAlign(22);
   // tex->SetTextSize(0.1);
   // tex->SetLineWidth(2);
   // tex->Draw();

   TCurlyLine *z0 = new TCurlyLine(30, 30, 110, 30, .04, .03);
   z0->SetLineWidth(2);
   z0->SetWavy();
   z0->Draw();

      tex = new TLatex(71,37.5,"#font[12]{Z^{0}}");
   tex->SetTextAlign(22);
   tex->SetTextSize(0.1);
   tex->SetLineWidth(2);
   tex->Draw();

   // line = new TLine(110,30,130,10);
   // line->SetLineWidth(3);
   // line->Draw();
   // line = new TLine(110,30,130,50);
   // line->SetLineWidth(3);
   // line->Draw();

   ar1 = new TArrow(110,30,130,10,0.03, "-<|-");
   ar1->SetLineWidth(3);
   ar1->SetFillColor(1);
   ar1->Draw();

   ar2 = new TArrow(110,30,130,50,0.03, "-|>-");
   ar2->SetLineWidth(3);
   ar2->SetFillColor(1);
   ar2->Draw();

      tex = new TLatex(135,6,"#font[12]{#bar{l}}");
   tex->SetTextAlign(22);
   tex->SetTextSize(0.1);
   tex->SetLineWidth(2);
   tex->Draw();

      tex = new TLatex(135,55,"#font[12]{l}");
   tex->SetTextAlign(22);
   tex->SetTextSize(0.1);
   tex->SetLineWidth(2);
   tex->Draw();

   c1->Modified();
   c1->cd();
}
